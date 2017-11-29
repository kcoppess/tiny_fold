#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <vector>
using namespace std;

string seqs[5] = {"AAAAA", "AAAAU", "GCCCC", "ACCCCGU", "AACCCCGUU"};
string circ[4] = {"AAAAA", "AAAAU", "ACCCUCCC", "AAACCCUUUCCC"};

vector<double> ener = {5.69, 6.0, 4.09, -7.09};

double tol = 1e-11;

bool almost_equal(double x, double y) {
    return abs(x-y) < tol;
}

TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[5] = {1,1.0000124791631115,1.0001857787828816,1.0021716139156089,1.0233027459871049};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition function for circular sequences", "[partition_circular]") {
    double part[4] = {1, 1, 1.0000023076971, 1.0003791934327};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition gradients for linear sequences", "[gradient_part_linear]") {
    vector< vector<double> > grad = {{0,0,0,0},
        {-2.10624588e-05,0,0,0},
        {0,0,-0.000313559325,0},
        {-0.00335171279,0,-0.00364420966,-0.003330650339},
        {-0.0746191942895,0,-0.039022635596,-0.074311205720273069}};
}