#include <iostream>
#include <string>
#include "rna.hh"
using namespace std;

int main() {
    string seq = "AUGC";
    vector<double> param;
    param.resize(4);
    for (int i=0; i < 4; i++) {
        param[i] = i;
    }
    RNA bob(seq, true, param);

    cout << bob.get_partition() << endl;
    bob.update_energy(param);
    return 0;
}
