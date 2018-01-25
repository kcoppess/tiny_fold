#include "rna.hh"
#include <iostream>
#include <string>
#include <valarray>
#include <vector>
#include <cmath>
#include <ctime>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

typedef std::vector< std::vector< std::valarray<double> > > tensor;
typedef std::vector< std::valarray<double> > matrix;
typedef std::valarray<double> vect;

/* Global variables */
const double R = 0.0019872; // kcal/K/mol universal gas constant
const double T = 298.15; // K temperature (room temp)
const double invRT = 1./(R*T);

const double g_loop = 1.; // kcal/mol
const double exp_neg_gloop_over_RT = exp(-invRT*g_loop);

/* Public Functions */

RNA::RNA(std::string seq, bool type, vect ener, bool wantBPP, bool grad) { 
    isCircular = type;
    calcBPP = wantBPP;
    calcGrad = grad;
    sequence = seq;
    nn = seq.length();
    energies = ener;
    
    g_base_pair.resize(nn, vect(nn));
    calc_gBasePair();
    
    partitionBound.resize(nn, vect(nn));
    partition.resize(nn, vect(nn));
    for (int ii = 1; ii < nn; ii++) {
        partition[ii][ii-1] = 1;
    }
    partitionS.resize(nn, vect(nn)); // storage matrix
    calc_partition();
    
    if (calcGrad) { 
        gradientBound.resize(nn, matrix(nn, vect(4)));
        gradientS.resize(nn, matrix(nn, vect(4))); // storage matrix
        gradient.resize(nn, matrix(nn, vect(4)));
        calc_gradient(); 
    }
    if (calcBPP) {
        bpp.resize(nn, vect(nn));
        bppS.resize(nn, vect(nn)); // storage matrix
        calc_bpp();
        if (calcGrad) { 
            bppGradientS.resize(nn, matrix(nn, vect(4))); // storage matrix
            bppGradient.resize(nn, matrix(nn, vect(4)));
            calc_bpp_gradient(); 
        }
    }
}

RNA::~RNA() {}

bool RNA::is_circular() { return isCircular; }

int RNA::get_length() { return nn; }

std::string RNA::get_sequence() { return sequence; }

vect RNA::get_energy() { return energies; }

double RNA::get_partition() { return partition[0][nn-1]; }

vect RNA::get_gradient() { return gradient[0][nn-1]; }

void RNA::update_energy(vect ener) {
    energies = ener;
    calc_gBasePair();
    calc_partition();
    if (calcGrad) { calc_gradient(); }
    if (calcBPP) {
        calc_bpp();
        if (calcGrad) { calc_bpp_gradient(); }
    }
}

matrix RNA::get_bpp_full() { return bpp; }

double RNA::get_bpp(int i, int j) { return bpp[i][j]; }

tensor RNA::get_bpp_gradient_full() { return bppGradient; }

vect RNA::get_bpp_gradient(int i, int j) { return bppGradient[i][j]; }

/* Private Functions */

void RNA::calc_gBasePair(){
    //g_base_pair.resize(nn, vect(nn));
    for (int ii=0; ii < nn; ii++) {
        for (int jj=0; jj < nn; jj++) {
            if ((sequence[ii] == 'A' && sequence[jj] == 'U') || (sequence[ii] == 'U' && sequence[jj] == 'A')) {
                g_base_pair[ii][jj] = energies[0];
            } else if ((sequence[ii] == 'U' && sequence[jj] == 'G') || (sequence[ii] == 'G' && sequence[jj] == 'U')) {
                g_base_pair[ii][jj] = energies[1];
            } else if ((sequence[ii] == 'G' && sequence[jj] == 'C') || (sequence[ii] == 'C' && sequence[jj] == 'G')) {
                g_base_pair[ii][jj] = energies[2];
            } else {
                g_base_pair[ii][jj] = 0.;
            }
        }
    }
}

int RNA::find_energy_index(int i, int j) {
    for (int mm = 0; mm < 4; mm++) {
        if (g_base_pair[i][j] == energies[mm]) {
            return mm;
        }
    }
    return -1;
}


// returns hairpin loop energy
double RNA::hairpin(double gHP) {
   return exp(-invRT * (gHP + g_loop));
}

double RNA::interior(double gBP, char loop) {
    if (loop == 's') {
        return exp(-invRT*(gBP + energies[3]));
    } else if (loop == 'l') {
        return exp(-invRT*(gBP + g_loop));
    } else{
        std::cout << loop << ": Invalid loop type" << std::endl;
        return 0;
    }
}

// calculates the partition function for every subsequence
void RNA::calc_partition() {
    //std::clock_t start = std::clock();
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));

    // entry i,j stores (bound) parition values for subsequence with ending bases i and j
    //partitionBound.resize(nn, vect(nn));
    //partition.resize(nn, vect(nn));
    //for (int ii = 1; ii < nn; ii++) {
    //    partition[ii][ii-1] = 1;
    //}
    //partitionS.resize(nn, vect(nn)); // storage matrix

    for (int ll = 1; ll < nn+1; ll++) { //iterating over all subsequence lengths
        for (int ii = 0; ii < nn-ll+1; ii++) { //iterating over all starting positions for subsequences
            int jj = ii + ll - 1; // ending position for subsequence
            
            // partitionBound recursion
            if (jj-ii > 3 && g_base_pair[ii][jj]) { // if possible hairpin: at least 4 positions apart and able to form a base pair
                if (isCircular) {
                    if ((ii + nn) - jj > 3) { // checking that base pair can form and bases are at least 4 positions apart on both sides
                        partitionBound[ii][jj] = hairpin(g_base_pair[ii][jj]);
                    }
                } else {
                    partitionBound[ii][jj] = hairpin(g_base_pair[ii][jj]);
                }
            } else { partitionBound[ii][jj] = 0; }
            for (int dd = ii+1; dd < jj-4; dd++) { // iterate over all possible rightmost pairs
                for (int ee = dd+4; ee < jj; ee++) { // i < d < e < j and d,e must be at least 4 positions apart
                    char interior_loop_type = ' ';
                    if (g_base_pair[ii][jj] && g_base_pair[dd][ee]) { // possible for both base pairs to form
                        if (ii+1 == dd && ee+1 == jj) { // if stacked
                            interior_loop_type = 's';
                        } else { // if loop
                            interior_loop_type = 'l';
                        }
                        if (isCircular) {
                            if ((dd + nn) - ee > 3) {
                                partitionBound[ii][jj] += partitionBound[dd][ee] * interior(g_base_pair[ii][jj], interior_loop_type);
                            }
                        } else {
                            partitionBound[ii][jj] += partitionBound[dd][ee] * interior(g_base_pair[ii][jj], interior_loop_type);
                        }
                    }
                }
            }

            partitionS[ii][jj] = 0.0;
            // partitionS recursion
            for (int dd = ii+4; dd < jj+1; dd++) { // iterate over all rightmost pairs with base i (beginning of subsequence)
                partitionS[ii][jj] += partitionBound[ii][dd];
            }

            // partition recursion
            partition[ii][jj] = 1.0;
            if (isCircular && ii == 0 && jj == nn-1) { // closing the chain for circular RNA
                for (int dd = 0; dd < nn-4; dd++) {
                    if (dd == 0) {
                        partition[ii][jj] += partitionS[0][jj]*exp_neg_gloop_over_RT;
                    } else {
                        if (partitionBound[0][dd-1] && partitionBound[dd][nn-1]) { // to account for stacked pair forming when chain is closed
                            partition[ii][jj] += (partition[0][dd-1] + partitionS[0][dd-1]*(exp_neg_gstack_gloop_over_RT - 1))*partitionS[dd][nn-1]*exp_neg_gloop_over_RT;
                        } else { // to account for interior loop forming when chain is closed
                            partition[ii][jj] += partition[ii][dd-1]*partitionS[dd][jj]*exp_neg_gloop_over_RT;
                        }
                    }
                }
            } else {
                for (int dd = ii; dd < jj-3; dd++) { // iterating over all possible rightmost pairs
                    if (dd == 0) { // to deal with issue of wrapping around in the last iteration
                        partition[ii][jj] += partitionS[dd][jj];
                    } else {
                        partition[ii][jj] += partition[ii][dd-1]*partitionS[dd][jj];
                    }
                }
            }
        }
    }
    //std::clock_t end = std::clock();
    //std::ofstream file("part.txt", std::ios_base::app);
    //file << 1e6 * (end-start)/CLOCKS_PER_SEC << " ";
    //file.close();
}

void RNA::calc_gradient() {
    //std::clock_t start = std::clock();
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));

    // entry i,j stores (bound) parition gradients for subsequence with ending bases i and j
    //gradientBound.resize(nn, matrix(nn, vect(4)));
    //gradientS.resize(nn, matrix(nn, vect(4))); // storage matrix
    //gradient.resize(nn, matrix(nn, vect(4)));

    for (int ll = 1; ll < nn+1; ll++) { //iterating over all subsequence lengths
        for (int ii = 0; ii < nn-ll+1; ii++) { //iterating over all starting positions for subsequences
            int jj = ii + ll - 1; // ending position for subsequence
            int ind = find_energy_index(ii, jj); // keeps track of which energy parameter
            
            gradientBound[ii][jj] = {0, 0, 0, 0};
            // partitionBound recursion
            if (jj-ii > 3 && g_base_pair[ii][jj]) { // if possible hairpin: at least 4 positions apart and able to form a base pair
                if (isCircular) {
                    if ((ii + nn) - jj > 3) { // checking that base pair can form and bases are at least 4 positions apart on both sides
                        gradientBound[ii][jj][ind] = -invRT*hairpin(g_base_pair[ii][jj]);
                    }
                } else {
                    gradientBound[ii][jj][ind] = -invRT*hairpin(g_base_pair[ii][jj]);
                }
            }
            for (int dd = ii+1; dd < jj-4; dd++) { // iterate over all possible rightmost pairs
                for (int ee = dd+4; ee < jj; ee++) { // i < d < e < j and d,e must be at least 4 positions apart
                    char interior_loop_type = ' ';
                    if (g_base_pair[ii][jj] && g_base_pair[dd][ee]) { // possible for both base pairs to form
                        if (ii+1 == dd && ee+1 == jj) { // if stacked
                            interior_loop_type = 's';
                        } else { // if loop
                            interior_loop_type = 'l';
                        }
                        if (isCircular) {
                            if ((dd + nn) - ee > 3) {
                                double inter = interior(g_base_pair[ii][jj], interior_loop_type);
                                gradientBound[ii][jj] += gradientBound[dd][ee] * inter;
                                gradientBound[ii][jj][ind] += -invRT * partitionBound[dd][ee] * inter;
                                if (interior_loop_type == 's') {
                                    gradientBound[ii][jj][3] += -invRT * inter * partitionBound[dd][ee];
                                }
                            }
                        } else {
                            double inter = interior(g_base_pair[ii][jj], interior_loop_type);
                            gradientBound[ii][jj] += gradientBound[dd][ee] * inter;
                            gradientBound[ii][jj][ind] += -invRT * partitionBound[dd][ee] * inter;
                            if (interior_loop_type == 's') {
                                gradientBound[ii][jj][3] += -invRT * inter * partitionBound[dd][ee];
                            }
                        }
                    }
                }
            }
            gradientS[ii][jj] = {0, 0, 0, 0};
            // partitionS recursion
            for (int dd = ii+4; dd < jj+1; dd++) { // iterate over all rightmost pairs with base i (beginning of subsequence)
                gradientS[ii][jj] += gradientBound[ii][dd];
            }
            gradient[ii][jj] = {0, 0, 0, 0};
            // partition recursion
            if (isCircular && ii == 0 && jj == nn-1) { // closing the chain for circular RNA
                for (int dd = 0; dd < nn-4; dd++) {
                    if (dd == 0) {
                        gradient[ii][jj] += gradientS[0][jj] * exp_neg_gloop_over_RT;
                    } else {
                        if (partitionBound[0][dd-1] && partitionBound[dd][nn-1]) { // to account for stacked pair forming when chain is closed
                            // chain rule
                            gradient[ii][jj] += gradient[0][dd-1] * partitionS[dd][nn-1] * exp_neg_gloop_over_RT;
                            gradient[ii][jj] += partition[0][dd-1] * gradientS[dd][nn-1] * exp_neg_gloop_over_RT;
                            gradient[ii][jj] += gradientS[0][dd-1] * (exp_neg_gstack_gloop_over_RT - 1) * partitionS[dd][nn-1] * exp_neg_gloop_over_RT;
                            gradient[ii][jj] += partitionS[0][dd-1] * (exp_neg_gstack_gloop_over_RT - 1) * gradientS[dd][nn-1] * exp_neg_gloop_over_RT;
                            gradient[ii][jj][3] += -invRT * partitionS[0][dd-1] * exp_neg_gstack_gloop_over_RT * partitionS[dd][nn-1] * exp_neg_gloop_over_RT;
                        } else { // to account for interior loop forming when chain is closed
                            gradient[ii][jj] += gradient[ii][dd-1] * partitionS[dd][jj] * exp_neg_gloop_over_RT;
                            gradient[ii][jj] += partition[ii][dd-1] * gradientS[dd][jj] * exp_neg_gloop_over_RT;
                        }
                    }
                }
            } else {
                for (int dd = ii; dd < jj-3; dd++) { // iterating over all possible rightmost pairs
                    if (dd == 0) { // to deal with issue of wrapping around in the last iteration
                        gradient[ii][jj] += gradientS[dd][jj];
                    } else {
                        gradient[ii][jj] += gradient[ii][dd-1] * partitionS[dd][jj] + partition[ii][dd-1] * gradientS[dd][jj];
                    }
                }
            }
        }
    }
    //std::clock_t end = std::clock();
    //std::ofstream file("grad.txt", std::ios_base::app);
    //file << 1e6 * (end-start)/CLOCKS_PER_SEC << " ";
    //file.close();
}

// XXX for circular RNA bpp calculations
void RNA::sum_left_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT) {
    for (int kk = 0; kk < ii-3; kk++) {
        for (int ll = kk+3; ll < ii; ll++) {
            if ((kk+nn == jj+1) && (ll == ii-1)) { //stacked pair
                bpp[ii][jj] += qbij_over_full * partitionBound[kk][ll] * exp_neg_gstack_over_RT; 
            } else { 
                bpp[ii][jj] += qbij_over_full * partitionBound[kk][ll] * exp_neg_gloop_over_RT;
            }
        }
    }
}

void RNA::gradient_sum_left_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT) {
    vect grad_full_part = gradient[0][nn-1];
    vect d_q_bound_ij = gradientBound[ii][jj];
    double term;
    double full_part = partition[0][nn-1];

    for (int kk = 0; kk < ii-3; kk++) {
        for (int ll = kk+3; ll < ii; ll++) {
            if ((kk+nn == jj+1) && (ll == ii-1)) { //stacked pair
                term = qbij_over_full * partitionBound[kk][ll] * exp_neg_gstack_over_RT; 
                bppGradient[ii][jj][3] += -invRT * term;
                bppGradient[ii][jj] += qbij_over_full * gradientBound[kk][ll] * exp_neg_gstack_over_RT;
                bppGradient[ii][jj] += d_q_bound_ij * partitionBound[kk][ll] * exp_neg_gstack_over_RT / full_part;
                bppGradient[ii][jj] += -term * grad_full_part / full_part;
            } else { 
                term = qbij_over_full * partitionBound[kk][ll] * exp_neg_gloop_over_RT; 
                bppGradient[ii][jj] += qbij_over_full * gradientBound[kk][ll] * exp_neg_gloop_over_RT;
                bppGradient[ii][jj] += d_q_bound_ij * partitionBound[kk][ll] * exp_neg_gloop_over_RT / full_part;
                bppGradient[ii][jj] += -term * grad_full_part / full_part;
            }
        }
    }
}

// XXX for circular RNA bpp calculations
void RNA::sum_right_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT) {
    //std::cout << ii << " " << jj << ": entered" << std::endl;
    for (int kk = jj+1; kk < nn-3; kk++) {
        for (int ll = kk+3; ll < nn; ll++) {
            //std::cout << kk << " " << ll << ": made it in" << std::endl;
            if ((kk == jj+1) && (ll+1 == ii+nn)) { //stacked pair
                bpp[ii][jj] += qbij_over_full * partitionBound[kk][ll] * exp_neg_gstack_over_RT; 
            } else { 
                bpp[ii][jj] += qbij_over_full * partitionBound[kk][ll] * exp_neg_gloop_over_RT;
            }
        }
    }
}

void RNA::gradient_sum_right_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT) {
    vect grad_full_part = gradient[0][nn-1];
    vect d_q_bound_ij = gradientBound[ii][jj];
    double term;
    double full_part = partition[0][nn-1];
    //std::cout << ii << " " << jj << ": entered" << std::endl;
    for (int kk = jj+1; kk < nn-3; kk++) {
        for (int ll = kk+3; ll < nn; ll++) {
            //std::cout << kk << " " << ll << ": made it in" << std::endl;
            if ((kk == jj+1) && (ll+1 == ii+nn)) { //stacked pair
                term = qbij_over_full * partitionBound[kk][ll] * exp_neg_gstack_over_RT; 
                bppGradient[ii][jj][3] += -invRT * term;
                bppGradient[ii][jj] += qbij_over_full * gradientBound[kk][ll] * exp_neg_gstack_over_RT;
                bppGradient[ii][jj] += d_q_bound_ij * partitionBound[kk][ll] * exp_neg_gstack_over_RT / full_part;
                bppGradient[ii][jj] += -term * grad_full_part / full_part;
            } else { 
                term = qbij_over_full * partitionBound[kk][ll] * exp_neg_gloop_over_RT; 
                bppGradient[ii][jj] += qbij_over_full * gradientBound[kk][ll] * exp_neg_gloop_over_RT;
                bppGradient[ii][jj] += d_q_bound_ij * partitionBound[kk][ll] * exp_neg_gloop_over_RT / full_part;
                bppGradient[ii][jj] += -term * grad_full_part / full_part;
            }
        }
    }
}

void RNA::sum_exterior_basepairs(int ii, int jj, double q_bound_ij, double exp_neg_gstack_over_RT) {
    for (int ll = 0; ll < ii; ll++) {
        bpp[ii][jj] += q_bound_ij * bppS[ll][jj];
    }
    if (partitionBound[ii-1][jj+1]) { // stacked pairs
        bpp[ii][jj] += -q_bound_ij * exp(-invRT * g_base_pair[ii-1][jj+1]) * bpp[ii-1][jj+1] * (exp_neg_gloop_over_RT - exp_neg_gstack_over_RT) / partitionBound[ii-1][jj+1];
    }
}

void RNA::gradient_sum_exterior_basepairs(int ii, int jj, double q_bound_ij, double exp_neg_gstack_over_RT) {
    vect d_q_bound_ij = gradientBound[ii][jj];
    for (int ll = 0; ll < ii; ll++) {
        bppGradient[ii][jj] += q_bound_ij * bppGradientS[ll][jj] + d_q_bound_ij * bppS[ll][jj];
    }
    if (partitionBound[ii-1][jj+1]) { // stacked pairs
        double bpp_ij = bpp[ii-1][jj+1];
        double q_denom = partitionBound[ii-1][jj+1];
        double exp_neg_invRT_gbp_ij = exp(-invRT * g_base_pair[ii-1][jj+1]);

        double term = -q_bound_ij * exp_neg_invRT_gbp_ij * bpp_ij * (exp_neg_gloop_over_RT - exp_neg_gstack_over_RT) / q_denom;
        bppGradient[ii][jj] += (-d_q_bound_ij * exp_neg_invRT_gbp_ij * bpp_ij * (exp_neg_gloop_over_RT - exp_neg_gstack_over_RT)) / q_denom;
        bppGradient[ii][jj] += -q_bound_ij * exp_neg_invRT_gbp_ij * bppGradient[ii-1][jj+1] * (exp_neg_gloop_over_RT - exp_neg_gstack_over_RT) / q_denom;
        bppGradient[ii][jj] += -term * gradientBound[ii-1][jj+1] / q_denom;
        
        int adj_energy_index = find_energy_index(ii-1, jj+1);
        bppGradient[ii][jj][adj_energy_index] += -invRT * term;
        bppGradient[ii][jj][3] += -invRT * q_bound_ij * exp_neg_invRT_gbp_ij * bpp_ij * exp_neg_gstack_over_RT / q_denom;
    }
}

void RNA::calc_bppS_entry(int ii, int jj) {
    bppS[ii][jj] = 0.0;
    for (int kk = jj+1; kk < nn; kk++) {
        double q_bound_ik = partitionBound[ii][kk];
        if (q_bound_ik) {
            bppS[ii][jj] += exp(-invRT * g_base_pair[ii][kk]) * exp_neg_gloop_over_RT * bpp[ii][kk] / q_bound_ik;
        }
    }
}

void RNA::calc_bppGradientS_entry(int ii, int jj) {
    bppGradientS[ii][jj] = {0, 0, 0, 0};
    for (int kk = jj+1; kk < nn; kk++) {
        double q_bound_ik = partitionBound[ii][kk];
        if (q_bound_ik) {
            int energy_index = find_energy_index(ii, kk);
            vect qbik_grad = gradientBound[ii][kk];
            double exp_neg_gbp_over_RT = exp(-invRT * g_base_pair[ii][kk]);
            double term = exp_neg_gbp_over_RT * exp_neg_gloop_over_RT * bpp[ii][kk] / q_bound_ik;
            
            bppGradientS[ii][jj][energy_index] += -invRT * term;
            bppGradientS[ii][jj] += exp_neg_gbp_over_RT * exp_neg_gloop_over_RT * bppGradient[ii][kk] / q_bound_ik;
            bppGradientS[ii][jj] += -term * qbik_grad / q_bound_ik;
        }
    }
}

void RNA::calc_bpp() {
    //std::clock_t start = std::clock();
    double exp_neg_gstack_over_RT = exp(-invRT*energies[3]);
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));
    double full_part = partition[0][nn-1];
    
    // base pair probability
    //bpp.resize(nn, vect(nn));
    //bppS.resize(nn, vect(nn)); // storage matrix
    
    // need to start with outside pairs and work way in
    for (int ii = 0; ii < nn; ii++) { // index for first base
        for (int jj = nn-1; jj > ii + 3; jj--) { // index for second base
            double q_bound_ij = partitionBound[ii][jj];
            double qbij_over_full = q_bound_ij / full_part;

            calc_bppS_entry(ii, jj);
            
            if (isCircular) {
                bpp[ii][jj] = qbij_over_full * exp(-invRT * g_base_pair[ii][jj]);
                sum_left_interior_loops(ii, jj, qbij_over_full, exp_neg_gstack_over_RT);
                sum_right_interior_loops(ii, jj, qbij_over_full, exp_neg_gstack_over_RT);
            } else {
                if (ii == 0 && jj == nn-1) {
                    bpp[ii][jj] = qbij_over_full;
                } else if (ii == 0) {
                    bpp[ii][jj] = qbij_over_full * partition[jj+1][nn-1];
                } else if (jj == nn-1) {
                    bpp[ii][jj] = partition[0][ii-1] * qbij_over_full;
                } else {
                    bpp[ii][jj] = partition[0][ii-1] * qbij_over_full * partition[jj+1][nn-1]; // if base-pair is not enclosed
                }
            }
            if (ii != 0 && jj != nn-1) {
                sum_exterior_basepairs(ii, jj, q_bound_ij, exp_neg_gstack_over_RT);
            }
        }
    }
    //std::clock_t end = std::clock();
    //std::ofstream file("bpp.txt", std::ios_base::app);
    //file << 1e6 * (end-start)/CLOCKS_PER_SEC << " ";
    //file.close();
}

// XXX NEED TO FIX circular sequence conditions
void RNA::calc_bpp_gradient() {
    //std::clock_t start = std::clock();
    double exp_neg_gstack_over_RT = exp(-invRT*energies[3]);
    double full_part = partition[0][nn-1];
    vect grad_full_part = gradient[0][nn-1];
    double term;
    
    //bppGradientS.resize(nn, matrix(nn, vect(4))); // storage matrix
    //bppGradient.resize(nn, matrix(nn, vect(4)));

    // need to start with outside pairs and work way in
    for (int ii = 0; ii < nn; ii++) { // index for first base
        for (int jj = nn-1; jj > ii+3; jj--) { // index for second base
            int energy_index = find_energy_index(ii, jj);
            
            double q_bound_ij = partitionBound[ii][jj];
            double qbij_over_full = q_bound_ij / full_part;
            vect d_q_bound_ij = gradientBound[ii][jj];
            double q_0i = partition[0][ii-1];
            double q_jN = 0;
            if (jj < nn-1) { q_jN = partition[jj+1][nn-1]; }

            calc_bppGradientS_entry(ii, jj);

            if (isCircular) {
                term = qbij_over_full * exp(-invRT * g_base_pair[ii][jj]);
                bppGradient[ii][jj] = d_q_bound_ij * exp(-invRT * g_base_pair[ii][jj]) / full_part;
                bppGradient[ii][jj] += -term * grad_full_part / full_part;
                bppGradient[ii][jj][energy_index] += -invRT * term;
                
                gradient_sum_left_interior_loops(ii, jj, qbij_over_full, exp_neg_gstack_over_RT);
                gradient_sum_right_interior_loops(ii, jj, qbij_over_full, exp_neg_gstack_over_RT);
            } else {
                if (ii == 0 && jj == nn-1) {
                    term = qbij_over_full;
                    bppGradient[ii][jj] = (d_q_bound_ij - term * grad_full_part) / full_part;
                } else if (ii == 0) {
                    term = qbij_over_full * q_jN;
                    bppGradient[ii][jj] = (d_q_bound_ij * q_jN) / full_part + qbij_over_full * gradient[jj+1][nn-1];
                    bppGradient[ii][jj] += -term * grad_full_part / full_part;
                } else if (jj == nn-1) {
                    term = q_0i * qbij_over_full;
                    bppGradient[ii][jj] = (gradient[0][ii-1] * qbij_over_full) + q_0i * d_q_bound_ij / full_part;
                    bppGradient[ii][jj] += -term * grad_full_part / full_part;
                } else {
                    term = q_0i * qbij_over_full * q_jN; // if base-pair is not enclosed
                    bppGradient[ii][jj] = (gradient[0][ii-1] * qbij_over_full * q_jN) + (q_0i * d_q_bound_ij * q_jN / full_part) + q_0i * qbij_over_full * gradient[jj+1][nn-1];
                    bppGradient[ii][jj] += -term * grad_full_part / full_part;
                }
            }
            if (ii != 0 && jj != nn-1) {
                gradient_sum_exterior_basepairs(ii, jj, q_bound_ij, exp_neg_gstack_over_RT);
            }
        }
    }
    //std::clock_t end = std::clock();
    //std::ofstream file("bpp_grad.txt", std::ios_base::app);
    //file << 1e6 * (end-start)/CLOCKS_PER_SEC << " ";
    //file.close();
}

PYBIND11_PLUGIN(tinyfold) {
    py::module m("tinyfold", "tinyfold plugin");
    py::class_<RNA>(m, "RNA")
        .def(py::init([] (std::string seq, bool type, vect ener, bool wantBPP, bool grad) {
                    return new RNA::RNA(seq, type, ener, wantBPP, grad);
                }))
        .def("is_circular", &RNA::is_circular, "indicates whether RNA is circular or not")
        .def("get_length", &RNA::get_length, "returns length of sequence (num of nucleotides)")
        .def("get_sequence", &RNA::get_sequence, "returns nucleotide sequence")
        .def("get_energy", &RNA::get_energy, "returns energy parameters")
        .def("get_partition", &RNA::get_partition, "returns value of partition function")
        .def("get_gradient", &RNA::get_gradient, "returns gradient of the partition function")
        .def("update_energy", &RNA::update_energy, "updates partition (and bpp and gradients) for new energy parameters")
        .def("get_bpp_full", &RNA::get_bpp_full, "returns full bpp matrix")
        .def("get_bpp", &RNA::get_bpp, "returns bpp for a particular pair of bases")
        .def("get_bpp_gradient_full", &RNA::get_bpp_gradient_full, "returns gradients for all basepair combinations")
        .def("get_bpp_gradient", &RNA::get_bpp_gradient, "returns gradient for particular pair of bases");
}
