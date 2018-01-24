#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <valarray>
using namespace std;

typedef vector< vector< valarray<double> > > tensor;
typedef vector< valarray<double> > matrix;
typedef valarray<double> vect;

string seqs[5] = {"AAAAA", "AAAAU", "GCCCC", "ACCCCGU", "AACCCCGUU"};
string circ[4] = {"AAAAA", "AAAAU", "ACCCUCCC", "AAACCCUUUCCC"};
string circ_1[4] = {"AAAAA", "AAUAA", "CUCCCACC", "ACCCUUUCCCAA"};
string circ_2[4] = {"AAAAA", "AUAAA", "CCACCCUC", "CCUUUCCCAAAC"};

vect ener = {5.69, 6.0, 4.09, -7.09};

double tol = 1e-11;

bool almost_equal(double x, double y) {
    return abs(x-y) < tol;
}

bool vect_almost_equal(vect x, vect y) {
    double diff2 = 0;
    for (int i = 0; i < 4; i++) {
        cout << x[i] << " " <<  y[i] << endl;
        diff2 += pow(x[i] - y[i],2);
    }
    return sqrt(diff2) < tol;
}
/*
TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[5] = {1,1.0000124791631115,1.0001857787828816,1.0021716139156089,1.0233027459871049};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener, false, false);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition function for circular sequences", "[partition_circular]") {
    double part[4] = {1, 1, 1.0000023076971, 1.0003791934327};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener, false, false);
        REQUIRE( almost_equal(bob.get_partition(), part[i]) );
    }
}

TEST_CASE("Partition gradients for linear sequences", "[gradient_part_linear]") {
    matrix grad = {{0,0,0,0},
        {-2.10624588e-05,0,0,0},
        {0,0,-0.000313559325,0},
        {-0.00335171279,0,-0.00364420966,-0.003330650339},
        {-0.0746191942895,0,-0.039022635596,-0.074311205720273069}};
    for (int i = 0; i < 5; i++) {
        RNA bob(seqs[i], false, ener, false, true);
        REQUIRE( vect_almost_equal(bob.get_gradient(), grad[i]) );
    } 
}

TEST_CASE("Partition gradients for circular sequences", "[gradient_part_circular]") {
    matrix grad = {{0,0,0,0},
        {0,0,0,0},
        {-3.8949547054e-06,0,0,0},
        {-0.0016844202181,0,0,-0.001044413431}};
    for (int i = 0; i < 4; i++) {
        RNA bob(circ[i], true, ener, false, true);
        REQUIRE( vect_almost_equal(bob.get_gradient(), grad[i]) );
    } 
}

TEST_CASE("Basepair probabilities for linear sequences", "[bpp_linear]") {
    matrix aaaaa(5, valarray<double>(5));
    matrix aaaau(5, valarray<double>(5));
    matrix gcccc(5, valarray<double>(5));
    matrix accccgu(7, valarray<double>(7));
    matrix aaccccguu(9, valarray<double>(9));
    
    aaaau[0][4] = 1.2479007383961551e-05;
    gcccc[0][4] = 0.00018574427553608532;
    accccgu[0][6] = 0.0019815320102396245;
    accccgu[1][5] = 0.0021544561056377051;
    aaccccguu[0][7] = 1.2197252009936157e-05;
    aaccccguu[0][8] = 0.020625537448534667;
    aaccccguu[1][7] = 0.022553953705591098;
    aaccccguu[1][8] = 1.2197252009936157e-05;
    aaccccguu[2][6] = 0.022593777679140784;

    matrix truth[] = {aaaaa, aaaau, gcccc, accccgu, aaccccguu};

    for (int i = 0; i < 5; i++) {
        RNA thor(seqs[i], false, ener, true, false);
        for (int j = 0; j < seqs[i].length(); j++) {
            REQUIRE( vect_almost_equal( thor.get_bpp_full()[j], truth[i][j]) );
        }
    }
}

TEST_CASE("Basepair probabilities for circular sequences", "[bpp_circular]") {
    matrix aaaaa(5, valarray<double>(5));
    matrix aaaau(5, valarray<double>(5));
    matrix acccuccc(8, valarray<double>(8));
    matrix aaacccuuuccc(12, valarray<double>(12));
    
    acccuccc[0][4] = 2.3076917718221965e-6;
    aaacccuuuccc[0][6] = 2.3068223653844651e-6;
    aaacccuuuccc[0][7] = 2.6810088901106141e-5;
    aaacccuuuccc[0][8] = 0.00028708534978963478;
    aaacccuuuccc[1][6] = 2.6810088901106141e-5;
    aaacccuuuccc[1][7] = 0.00031158850117650618;
    aaacccuuuccc[1][8] = 2.6810088901106141e-5;
    aaacccuuuccc[2][6] = 0.00028708534978963478;
    aaacccuuuccc[2][7] = 2.6810088901106141e-5;
    aaacccuuuccc[2][8] = 2.3068223653844651e-6;

    matrix truth[] = {aaaaa, aaaau, acccuccc, aaacccuuuccc};

    for (int i = 0; i < 4; i++) {
        RNA thor(circ[i], true, ener, true, false);
        for (int j = 0; j < circ[i].length(); j++) {
            REQUIRE( vect_almost_equal( thor.get_bpp_full()[j], truth[i][j]) );
        }
    }
}
*/
TEST_CASE("Gradient of basepair probabilities for linear sequences","[grad_bpp_linear]") {
    tensor aaaaa(5, matrix(5, vect(4)));
    tensor aaaau(5, matrix(5, vect(4)));
    tensor gcccc(5, matrix(5, vect(4)));
    tensor accccgu(7, matrix(7, vect(4)));
    tensor aaccccguu(9, matrix(9, vect(4)));
    
    //AU
    aaaau[0][4][0] = -2.1061933146135372e-05;
    accccgu[0][6][0] = -0.0033378228091143717;
    accccgu[1][5][0] = -0.0033162276547285432;
    aaccccguu[0][7][0] = -1.9697223269939332e-05;
    aaccccguu[0][8][0] = -0.068099385089383582;
    aaccccguu[1][7][0] = -0.07121356568856102;
    aaccccguu[1][8][0] = -1.9697223269939332e-05;
    aaccccguu[2][6][0] = -0.070752825849043341;

    //GC
    gcccc[0][4][2] = -0.00031344285229964632;
    accccgu[0][6][2] = -0.0033162276547285428;
    accccgu[1][5][2] = -0.0036284787194924687;
    aaccccguu[0][7][2] = 4.6130625989703884e-07;
    aaccccguu[0][8][2] = -0.033786243034938615;
    aaccccguu[1][7][2] = -0.036967505426624511;
    aaccccguu[1][8][2] = 4.6130625989703884e-07;
    aaccccguu[2][6][2] = -0.037272417173290488;

    //stack
    accccgu[0][6][3] = -0.0033168476362080901;
    accccgu[1][5][3] = -0.0033162729345412763;
    aaccccguu[0][7][3] = 8.8575204834229921e-07;
    aaccccguu[0][8][3] = -0.067866373951630518;
    aaccccguu[1][7][3] = -0.070981138779790301;
    aaccccguu[1][8][3] = 8.8575204834229921e-07;
    aaccccguu[2][6][3] = -0.070759614540283317;

    tensor truth[] = {aaaaa, aaaau, gcccc, accccgu, aaccccguu};

    for (int i = 0; i < 5; i++) {
        RNA thor(seqs[i], false, ener, true, true);
        for (int j = 0; j < seqs[i].length(); j++) {
            for (int k = 0; k < seqs[i].length(); k++) {
                REQUIRE( vect_almost_equal( thor.get_bpp_gradient_full()[j][k], truth[i][j][k]) );
            }
        }
    }

}

TEST_CASE("Gradient of basepair probabilities for circular sequences","[grad_bpp_circular]") {
    tensor aaaaa(5, matrix(5, vect(4)));
    tensor aaaau(5, matrix(5, vect(4)));
    tensor acccuccc(8, matrix(8, vect(4)));
    tensor aaacccuuuccc(12, matrix(12, vect(4)));
    
    //AU
    acccuccc[0][4][0] = -3.8949367286891062e-06;
    aaacccuuuccc[0][6][0] = -3.8895941385932963e-06;
    aaacccuuuccc[0][7][0] = -8.6562062499709132e-05;
    aaacccuuuccc[0][8][0] = -0.0014040094410586953;
    aaacccuuuccc[1][7][0] = -0.0014866815209148882;
    aaacccuuuccc[2][8][0] = -3.8895941385932963e-06;
    aaacccuuuccc[1][6][0] = -8.6562062499709132e-05;
    aaacccuuuccc[1][8][0] = -8.6562062499709132e-05;
    aaacccuuuccc[2][7][0] = -8.6562062499709132e-05;
    aaacccuuuccc[2][6][3] = -4.1328824466355418e-05;

    
    //stack
    aaacccuuuccc[0][6][3] = 2.4083630249027296e-09;
    aaacccuuuccc[0][7][3] = -4.1328824466355418e-05;
    aaacccuuuccc[0][8][3] = -0.0009196473798212065;
    aaacccuuuccc[1][7][3] = -0.00096097861277080429;
    aaacccuuuccc[2][8][3] = 2.4083630249027296e-09;
    aaacccuuuccc[1][6][3] = -4.1328824466355418e-05;
    aaacccuuuccc[1][8][3] = -4.1328824466355418e-05;
    aaacccuuuccc[2][7][3] = -4.1328824466355418e-05;
    aaacccuuuccc[2][6][3] = -0.0009196473798212065;

    tensor truth[] = {aaaaa, aaaau, acccuccc, aaacccuuuccc};

    for (int i = 0; i < 4; i++) {
        RNA thor(circ[i], true, ener, true, true);
        for (int j = 0; j < circ[i].length(); j++) {
            for (int k = 0; k < circ[i].length(); k++) {
                //cout << j << " " << k << endl;
                REQUIRE( vect_almost_equal( thor.get_bpp_gradient_full()[j][k], truth[i][j][k]) );
            }
        }
    }

}

