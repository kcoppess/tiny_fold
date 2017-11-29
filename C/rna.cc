#include "rna.hh"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

typedef vector< vector<double> > matrix;

/* Global variables */
const double R = 0.0019872; // kcal/K/mol universal gas constant
const double T = 298.15; // K temperature (room temp)
const double invRT = 1./(R*T);

const double g_loop = 1.; // kcal/mol
const double exp_neg_gloop_over_RT = exp(-invRT*g_loop);

/* Public Functions */

RNA::RNA(string seq, bool type, vector<double> ener) { 
    isCircular = type;
    sequence = seq;
    nn = seq.length();
    energies = ener;
    calc_partition();
}

RNA::~RNA() {}

bool RNA::is_circular() { return isCircular; }

int RNA::get_length() { return nn; }

string RNA::get_sequence() { return sequence; }

double RNA::get_energy() { return energies[3]; }

double RNA::get_partition() { return partition[0][nn-1]; }

void RNA::update_energy(vector<double> ener) {
    energies = ener;
    calc_partition();
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
        cout << loop << ": Invalid loop type" << endl;
        return 0;
    }
}

// calculates the partition function for every subsequence
void RNA::calc_partition() {
    double exp_neg_gstack_gloop_over_RT = exp(-invRT*(energies[3] - g_loop));

    // stores energies for each possible 
    matrix g_base_pair(nn, vector<double>(nn));
    calc_gBasePair(g_base_pair);

    // entry i,j stores (bound) parition values for subsequence with ending bases i and j
    partitionBound.resize(nn, vector<double>(nn));
    partition.resize(nn, vector<double>(nn));
    for (int ii = 0; ii < nn; ii++) {
        partition[ii][ii-1] = 1;
    }
    matrix partitionS(nn, vector<double>(nn)); // storage matrix

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

