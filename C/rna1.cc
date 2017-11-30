#include "rna.hh"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

typedef std::vector< std::vector<double> > matrix;

/* Global variables */
const double R = 0.0019872; // kcal/K/mol universal gas constant
const double T = 298.15; // K temperature (room temp)
const double invRT = 1./(R*T);

const double g_loop = 1.; // kcal/mol
const double exp_neg_gloop_over_RT = exp(-invRT*g_loop);

/* Public Functions */

RNA::RNA(std::string seq, bool type, std::vector<double> ener) { 
    isCircular = type;
    sequence = seq;
    nn = seq.length();
    energies = ener;
    calc_partition();
    calc_gradient();
}

RNA::~RNA() {}

bool RNA::is_circular() { return isCircular; }

int RNA::get_length() { return nn; }

std::string RNA::get_sequence() { return sequence; }

double RNA::get_energy() { return energies[3]; }

double RNA::get_partition() { return partition[0][nn-1]; }

std::vector<double> RNA::get_gradient() { return gradient[0][nn-1]; }

void RNA::update_energy(std::vector<double> ener) {
    energies = ener;
    calc_partition();
    calc_gradient();
}

/* Private Functions */

void RNA::calc_gBasePair(matrix& gBP){
    for (int ii=0; ii < nn; ii++) {
        for (int jj=0; jj < nn; jj++) {
            if ((sequence[ii] == 'A' && sequence[jj] == 'U') || (sequence[ii] == 'U' && sequence[jj] == 'A')) {
                gBP[ii][jj] = energies[0];
            } else if ((sequence[ii] == 'U' && sequence[jj] == 'G') || (sequence[ii] == 'G' && sequence[jj] == 'U')) {
                gBP[ii][jj] = energies[1];
            } else if ((sequence[ii] == 'G' && sequence[jj] == 'C') || (sequence[ii] == 'C' && sequence[jj] == 'G')) {
                gBP[ii][jj] = energies[2];
            } else {
                gBP[ii][jj] = 0.;
            }
        }
    }
    return;
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
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));

    // stores energies for each possible 
    matrix g_base_pair(nn, std::vector<double>(nn));
    calc_gBasePair(g_base_pair);

    // entry i,j stores (bound) parition values for subsequence with ending bases i and j
    partitionBound.resize(nn, std::vector<double>(nn));
    partition.resize(nn, std::vector<double>(nn));
    for (int ii = 0; ii < nn; ii++) {
        partition[ii][ii-1] = 1;
    }
    matrix partitionS(nn, std::vector<double>(nn)); // storage matrix

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
                                partitionBound[ii][jj] += partitionBound[dd][ee] * interior(g_base_pair[ii][jj], interior_loop_type);
                            }
                        } else {
                            partitionBound[ii][jj] += partitionBound[dd][ee] * interior(g_base_pair[ii][jj], interior_loop_type);
                        }
                    }
                }
            }
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
}

void RNA::calc_gradient() {
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));
    vector<double> grad(4);

    // stores energies for each possible 
    matrix g_base_pair(nn, std::vector<double>(nn));
    calc_gBasePair(g_base_pair);

    // entry i,j stores (bound) parition gradients for subsequence with ending bases i and j
    gradientBound.resize(nn, matrix(nn, std::vector<double>(4)));
    std::vector<matrix> gradientS(nn, matrix(nn, std::vector<double>(4))); // storage matrix
    gradient.resize(nn, matrix(nn, std::vector<double>(4)));

    for (int ll = 1; ll < nn+1; ll++) { //iterating over all subsequence lengths
        for (int ii = 0; ii < nn-ll+1; ii++) { //iterating over all starting positions for subsequences
            int jj = ii + ll - 1; // ending position for subsequence
            int ind = std::find(energies.begin(), energies.end(), g_base_pair[ii][jj]) - energies.begin();
            
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
                                grad = gradientBound[dd][ee];
                                std::transform(grad.begin(), grad.end(), grad.begin(), 
                                        std::bind1st(std::multiplies<double>, inter)); // gradientBound[dd][ee] * inter
                                // gradientBound[ii][jj] += gradientBound[dd][ee] * inter
                                std::transform(gradientBound[ii][jj].begin(), gradientBound[ii][jj].end(), grad.begin(), 
                                        gradientBound[ii][jj].begin(), std::plus<double>());
                                gradientBound[ii][jj][ind] += -invRT * partitionBound[dd][ee] * inter;
                            }
                        } else {
                            double inter = interior(g_base_pair[ii][jj], interior_loop_type);
                            grad = gradientBound[dd][ee];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, inter)); // gradientBound[dd][ee] * inter
                            // gradientBound[ii][jj] += gradientBound[dd][ee] * inter
                            std::transform(gradientBound[ii][jj].begin(), gradientBound[ii][jj].end(), grad.begin(), 
                                    gradientBound[ii][jj].begin(), std::plus<double>());
                            gradientBound[ii][jj][ind] += -invRT * partitionBound[dd][ee] * inter;
                        }
                    }
                }
            }
            // partitionS recursion
            for (int dd = ii+4; dd < jj+1; dd++) { // iterate over all rightmost pairs with base i (beginning of subsequence)
                // gradientS[ii][jj] += gradientBound[ii][dd]
                std::transform(gradientS[ii][jj].begin(), gradientS[ii][jj].end(), gradientBound[ii][dd].begin(), 
                        gradientS[ii][jj].begin(), std::plus<double>());
            }
            // partition recursion
            if (isCircular && ii == 0 && jj == nn-1) { // closing the chain for circular RNA
                for (int dd = 0; dd < nn-4; dd++) {
                    if (dd == 0) {
                        grad = gradientS[0][jj];
                        std::transform(grad.begin(), grad.end(), grad.begin(), 
                                std::bind1st(std::multiplies<double>, exp_neg_gloop_over_RT));
                        // gradient[ii][jj] += gradientS[0][jj] * exp_neg_gloop_over_RT
                        std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                gradient[ii][jj].begin(), std::plus<double>());
                    } else {
                        if (partitionBound[0][dd-1] && partitionBound[dd][nn-1]) { // to account for stacked pair forming when chain is closed
                            // gradient[ii][jj] += gradient[0][dd-1] * partitionS[dd][nn-1] * exp_neg_gloop_over_RT 
                            grad = gradient[0][dd-1];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, partitionS[dd][nn-1] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                            // gradient[ii][jj] += partition[0][dd-1] * gradientS[dd][nn-1] * exp_neg_gloop_over_RT
                            grad = gradientS[dd][nn-1];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, partition[0][dd-1] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                            // gradient[ii][jj] += gradientS[0][dd-1] * (exp_neg_gstack_gloop_over_RT - 1) * partitionS[dd][nn-1] * exp_neg_gloop_over_RT
                            double x = exp_neg_gstack_gloop_over_RT - 1;
                            grad = gradientS[0][dd-1];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, x * partitionS[dd][nn-1] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                            // gradient[ii][jj] += partitionS[0][dd-1] * (exp_neg_gstack_gloop_over_RT - 1) * gradientS[dd][nn-1] * exp_neg_gloop_over_RT
                            grad = gradientS[dd][nn-1];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, x * partitionS[0][dd-1] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                            
                            gradient[ii][jj][3] += -invRT * partitionS[0][dd-1] * exp_neg_gstack_gloop_over_RT * partitionS[dd][nn-1] * exp_neg_gloop_over_RT
                        } else { // to account for interior loop forming when chain is closed
                            // gradient[ii][jj] += gradient[ii][dd-1] * partitionS[dd][jj] * exp_neg_gloop_over_RT
                            grad = gradient[ii][dd-1];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, partitionS[dd][jj] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                            // gradient[ii][jj] += partition[ii][dd-1] * gradientS[dd][jj] * exp_neg_gloop_over_RT
                            grad = gradientS[dd][jj];
                            std::transform(grad.begin(), grad.end(), grad.begin(), 
                                    std::bind1st(std::multiplies<double>, partition[ii][dd-1] * exp_neg_gloop_over_RT));
                            std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                    gradient[ii][jj].begin(), std::plus<double>());
                        }
                    }
                }
            } else {
                for (int dd = ii; dd < jj-3; dd++) { // iterating over all possible rightmost pairs
                    if (dd == 0) { // to deal with issue of wrapping around in the last iteration
                        // gradient[ii][jj] += gradientS[dd][jj]
                        std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), gradientS[dd][jj].begin(), 
                                gradient[ii][jj].begin(), std::plus<double>());
                    } else {
                        // XXX finished editing here XXX XXX XXX XXX XXX
                        // gradient[ii][jj] += gradient
                        grad = gradientS[dd][jj];
                        std::transform(grad.begin(), grad.end(), grad.begin(), 
                                std::bind1st(std::multiplies<double>, partition[ii][dd-1] * exp_neg_gloop_over_RT));
                        std::transform(gradient[ii][jj].begin(), gradient[ii][jj].end(), grad.begin(), 
                                gradient[ii][jj].begin(), std::plus<double>());
                        partition[ii][jj] += partition[ii][dd-1]*partitionS[dd][jj];
                    }
                }
            }
        }
    }
}
