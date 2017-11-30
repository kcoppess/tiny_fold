#include <iostream>
#include <string>
#include <vector>
#include "rna.hh"
using namespace std;

typedef vector<vector<double> > matrix;

int main() {
    string seq = "AAAUUUUAAA";
    vector<double> ener = {5.69, 6.0, 4.09, -7.09};
    
    RNA bob(seq, false, ener);
    cout << bob.get_gradient()[0] << bob.get_gradient()[1] << bob.get_gradient()[2] << bob.get_gradient()[3] << endl;
    return 0;
}
