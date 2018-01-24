#include <ctime>
#include "rna.hh"
#include <string>
#include <iostream>
#include <valarray>

//XXX gotta comment out calc_gradient() in constructor to get correct timing of partition calculation

static const char bases[] = "AUGC";
const int L = 500;

std::valarray<double> ener = {5.69, 6.0, 4.09, -7.09};

int main() {
    std::srand(time(0));
    std::string sequence;
    for (int z = 0; z < L; z++) {
        sequence += bases[ rand() % 4];
    }
    std::cout << sequence << std::endl;
    //std::clock_t start = std::clock();
    RNA thor(sequence, false, ener, false, false);
    //std::clock_t end = std::clock();
    //std::cout << (end-start)/CLOCKS_PER_SEC << std::endl;

    return 0;
}
