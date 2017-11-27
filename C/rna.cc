#include "rna.hh"
#include <iostream>
#include <string>
#include <vector>
using namespace std;

typedef vector<vector<double>> matrix;

RNA::RNA(string seq, bool type, vector<double> ener) { 
    g_loop = 1.0;
    isCircular = type;
    sequence = seq;
    nn = seq.length();
    energies = ener;
    calc_partitionBound();
    calc_partition();
}

RNA::~RNA() {}

bool RNA::is_circular() { return isCircular; }

int RNA::get_length() { return nn; }

string RNA::get_sequence() { return sequence; }

double RNA::get_energy() { return energies[3]; }

double RNA::get_partition() { return partition[nn-1][nn-1]; }

void RNA::update_energy(vector<double> ener) {
    energies = ener;
    calc_partitionBound();
    calc_partition();
}

void RNA::calc_partitionBound() {
    partitionBound.resize(nn);
    for (int ii = 0; ii < nn; ii++) {
        partitionBound[ii].resize(nn);
    }
}

void RNA::calc_partition() {
    partition.resize(nn);
    for (int ii = 0; ii < nn; ii++) {
        partition[ii].resize(nn);
    }
}

void RNA::calc_gBasePair(matrix& gBP){
    return;
}
