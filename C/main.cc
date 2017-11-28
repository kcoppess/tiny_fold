#include <iostream>
#include <string>
#include <vector>
#include "rna.hh"
using namespace std;

typedef vector<vector<double> > matrix;

int main() {
    string seq = "ACCCUCCC";
    vector<double> ener = {5.69, 6.0, 4.09, -7.09};
    
    RNA bob(seq, true, ener);
    cout << bob.get_partition() << endl;
    return 0;
}
