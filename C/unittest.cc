#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <vector>
using namespace std;

string seqs[5] = {"AAAAA", "AAAAU", "GCCCC", "ACCCCGU", "AACCCCGUU"};
string circ[4] = {"AAAAA", "AAAAU", "ACCCUCCC", "AAACCCUUUCCC"};

vector<double> ener = {5.69, 6.0, 4.09, -7.09};

TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[5] = {1,1.0000124791631115,1.0001857787828816,1.0021716139156089,1.0233027459871049};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener);
        REQUIRE( bob.get_partition()  == part[i] );
    }
}

TEST_CASE("Partition function for circular sequences", "[partition_circular]") {
    double part[4] = {1, 1, 1.0000023076971, 1.0003791934327};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener);
        cout << bob.get_partition() - part[i] << ' ' << circ[i] << endl;
        REQUIRE( bob.get_partition() == part[i] );
    }
}
