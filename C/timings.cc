#include <ctime>
#include <fstream>
#include "rna.hh"
#include <string>
#include <iostream>
#include <valarray>

static const char bases[] = "AUGC";
const int L = 50;
const int num = 20;

std::valarray<double> ener = {5.69, 6.0, 4.09, -7.09};

int main() {
    std::ofstream file;
    std::srand(time(0));
    
    for (int l = 5; l < L+1; l++) {
        for (int i = 0; i < 20; i++) {
            std::string sequence;
            for (int z = 0; z < l; z++) {
                sequence += bases[ rand() % 4];
            }
            RNA thor(sequence, false, ener, true);
        }
        file.open("part.txt", std::ios_base::app);
        file << "\n";
        file.close();
        file.open("grad.txt", std::ios_base::app);
        file << "\n";
        file.close();
        file.open("bpp.txt", std::ios_base::app);
        file << "\n";
        file.close();
        file.open("bpp_grad.txt", std::ios_base::app);
        file << "\n";
        file.close();
    }

    //file.open("example.txt", std::ios_base::app);
    //file << "hello ";
    //file.close();
    return 0;
}
