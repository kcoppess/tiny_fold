#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include "rna.hh"
using namespace std;

typedef std::vector< std::vector< std::valarray<double> > > tensor;
typedef std::vector< std::valarray<double> > matrix;
typedef std::valarray<double> vect;

int main() {
    string seq = "AAAAU";
    double energy[] = {5.69, 6.0, 4.09, -7.09};
    vect ener(energy, 1);
    
    RNA bob(seq, false, ener);
    cout << bob.get_energy()[0] << endl;
    //cout << bob.get_gradient()[0] << bob.get_gradient()[1] << bob.get_gradient()[2] << bob.get_gradient()[3] << endl;
    return 0;
}
