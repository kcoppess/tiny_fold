#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <valarray>
using namespace std;

string seqs[5] = {"AAAAA", "AAAAU", "GCCCC", "ACCCCGU", "AACCCCGUU"};
string circ[4] = {"AAAAA", "AAAAU", "ACCCUCCC", "AAACCCUUUCCC"};

valarray<double> ener = {5.69, 6.0, 4.09, -7.09};

double tol = 1e-11;

bool almost_equal(double x, double y) {
    return abs(x-y) < tol;
}

bool vect_almost_equal(valarray<double> x, valarray<double> y) {
    double diff2 = 0;
    for (int i = 0; i < 4; i++) {
        diff2 += pow(x[i] - y[i],2);
    }
    return sqrt(diff2) < tol;
}

TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[5] = {1,1.0000124791631115,1.0001857787828816,1.0021716139156089,1.0233027459871049};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener, false);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition function for circular sequences", "[partition_circular]") {
    double part[4] = {1, 1, 1.0000023076971, 1.0003791934327};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener, false);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition gradients for linear sequences", "[gradient_part_linear]") {
    vector< valarray<double> > grad = {{0,0,0,0},
        {-2.10624588e-05,0,0,0},
        {0,0,-0.000313559325,0},
        {-0.00335171279,0,-0.00364420966,-0.003330650339},
        {-0.0746191942895,0,-0.039022635596,-0.074311205720273069}};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener, false);
        REQUIRE( vect_almost_equal(bob.get_gradient(), grad[i]) );
    } 
}

TEST_CASE("Partition gradients for circular sequences", "[gradient_part_circular]") {
    vector< valarray<double> > grad = {{0,0,0,0},
        {0,0,0,0},
        {-3.8949547054e-06,0,0,0},
        {-0.0016844202181,0,0,-0.001044413431}};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener, false);
        REQUIRE( vect_almost_equal(bob.get_gradient(), grad[i]) );
    } 
}

TEST_CASE("Basepair probabilities for linear sequences", "[bpp_linear]") {
    vector< valarray<double> > aaaaa(5, valarray<double>(5));
    vector< valarray<double> > aaaau(5, valarray<double>(5));
    vector< valarray<double> > gcccc(5, valarray<double>(5));
    vector< valarray<double> > accccgu(7, valarray<double>(7));
    vector< valarray<double> > aaccccguu(9, valarray<double>(9));
    
    aaaau[0][4] = 1.2479007383961551e-05;
    gcccc[0][4] = 0.00018574427553608532;
    accccgu[0][6] = 0.0019815320102396245;
    accccgu[1][5] = 0.0021544561056377051;
    aaccccguu[0][7] = 1.2197252009936157e-05;
    aaccccguu[0][8] = 0.020625537448534667;
    aaccccguu[1][7] = 0.022553953705591098;
    aaccccguu[1][8] = 1.2197252009936157e-05;
    aaccccguu[2][6] = 0.022593777679140784;

    vector< valarray<double> > truth[] = {aaaaa, aaaau, gcccc, accccgu, aaccccguu};

    for (int i = 0; i < 5; i++) {
        RNA thor(seqs[i], false, ener, true);
        for (int j = 0; j < seqs[i].length(); j++) {
            REQUIRE( vect_almost_equal( thor.get_bpp_full()[j], truth[i][j]) );
        }
    }
}
/*
def test_bpp_derivs_AU():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = -2.1061933146135372e-05
    gcccc = np.zeros((5,5))
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033378228091143717
    accccgu[1,5] = -0.0033162276547285432
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = -1.9697223269939332e-05
    aaccccguu[0,8] = -0.068099385089383582
    aaccccguu[1,7] = -0.07121356568856102
    aaccccguu[1,8] = -1.9697223269939332e-05
    aaccccguu[2,6] = -0.070752825849043341
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_AU(seq), val, tol))

def test_bpp_deriv_GC():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    gcccc = np.zeros((5,5))
    gcccc[0,4] = -0.00031344285229964632
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033162276547285428
    accccgu[1,5] = -0.0036284787194924687
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 4.6130625989703884e-07
    aaccccguu[0,8] = -0.033786243034938615
    aaccccguu[1,7] = -0.036967505426624511
    aaccccguu[1,8] = 4.6130625989703884e-07
    aaccccguu[2,6] = -0.037272417173290488
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_GC(seq), val, tol))

def test_bpp_deriv_stack():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    gcccc = np.zeros((5,5))
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033168476362080901
    accccgu[1,5] = -0.0033162729345412763
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 8.8575204834229921e-07
    aaccccguu[0,8] = -0.067866373951630518
    aaccccguu[1,7] = -0.070981138779790301
    aaccccguu[1,8] = 8.8575204834229921e-07
    aaccccguu[2,6] = -0.070759614540283317
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_stack(seq), val, tol))
*/
