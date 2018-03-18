#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <valarray>
using namespace std;

typedef vector< vector< valarray<double> > > tensor;
typedef vector< valarray<double> > matrix;
typedef valarray<double> vect;

string seqs[20] = {"AAAAA", "AUAAA", "AAAAU", "GAAAC", "GCCCC", "AACCCUU", "AUCCCAU", "UACCCUA", "CUAAAAG",
                    "CAAAAUG", "GUAAAAC", "GAAAAUC", "CGAAACG", "GGAAACC", "GCAAAGC", "AAUUUUU", 
                    "AAUUUUUUU", "UAUUUUUUA", "GGCCCCCCC", "CGCCCCCCG"};

vect ener = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2};

double tol = 1e-11;

bool almost_equal(double x, double y) {
    cout << x << " " << y << endl;
    return abs(x-y) < tol;
}

bool vect_almost_equal(vect x, vect y) {
    double diff2 = 0;
    for (int i = 0; i < 12; i++) {
        //cout << x[i] << " " <<  y[i] << endl;
        diff2 += pow(x[i] - y[i],2);
    }
    return sqrt(diff2) < tol;
}

TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[20] = {1, 1, 1.003811342961177, 1.028885907702594, 1.028885907702594, 1.015670157043195, 1.0079814993687743, 
                        1.007925773435268, 1.0346375796305893, 1.0343362346167273, 1.0340816903355927, 1.0338666783457544, 
                        1.0652583475428457, 1.1218674584946227, 1.0631135137836523, 1.0194815000043722, 1.0356200212515567,
                        1.0270115409114888, 1.2814478393674946, 1.2113566773834208};
    for (int i = 0; i < 19; i++) {
        RNA loki(seqs[i], false, ener, false, false);
        cout << seqs[i] << endl;
        REQUIRE( almost_equal(loki.get_partition(), part[i]) );
    }
}

/*
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

*/
