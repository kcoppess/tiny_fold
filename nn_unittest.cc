#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "rna.hh"
#include <string>
#include <valarray>
using namespace std;

typedef vector< vector< valarray<double> > > tensor;
typedef vector< valarray<double> > matrix;
typedef valarray<double> vect;

// testing for T = 298.15 K and g_loop = 1. kcal/mol
string seqs[21] = {"AAAAA", "AUAAA", "AAAAU", "GAAAC", "GCCCC", "AACCCUU", "AUCCCAU", "UACCCUA", "CUAAAAG",
                    "CAAAAUG", "GUAAAAC", "GAAAAUC", "CGAAACG", "GGAAACC", "GCAAAGC", "AAUUUUU", 
                    "AAUUUUUUU", "UAUUUUUUA", "GGCCCCCCC", "CGCCCCCCG", "AAACCCUUU"};

vect ener = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2};

double tol = 1e-11;

bool almost_equal(double x, double y) {
    cout << abs(x - y) << endl;
    return abs(x-y) < tol;
}

bool vect_almost_equal(vect x, vect y) {
    double diff2 = 0;
    for (int i = 0; i < y.size(); i++) {
        diff2 += pow(x[i] - y[i],2);
        cout << "  " << x[i] << " " << y[i] << endl;
    }
    //cout << diff2 << endl;
    return sqrt(diff2) < tol;
}

TEST_CASE("Partition function for linear sequences", "[partition_linear]") {
    double part[21] = {1, 1, 1.003811342961177, 1.028885907702594, 1.028885907702594, 1.015670157043195, 1.0079814993687743, 
                        1.007925773435268, 1.0346375796305893, 1.0343362346167273, 1.0340816903355927, 1.0338666783457544, 
                        1.0652583475428457, 1.1218674584946227, 1.0631135137836523, 1.0194815000043722, 1.0356200212515567,
                        1.0270115409114888, 1.2814478393674946, 1.2113566773834208, 1.0364471989019666}; //1.0411293813403846
    for (int i = 0; i < 21; i++) {
        cout << seqs[i] << endl;
        RNA loki(seqs[i], false, ener, false, false);
        REQUIRE( almost_equal(loki.get_partition(), part[i]) );
    }
}

/*
TEST_CASE("Basepair probabilities for linear sequences", "[bpp_linear]") {
    matrix aaaaa(5, valarray<double>(5));
    matrix auaaa(5, valarray<double>(5));
    matrix aaaau(5, valarray<double>(5));
    matrix gaaac(5, valarray<double>(5));
    matrix gcccc(5, valarray<double>(5));
    matrix aacccuu(7, valarray<double>(7));
    matrix aucccau(7, valarray<double>(7));
    matrix uacccua(7, valarray<double>(7));
    matrix cuaaaag(7, valarray<double>(7));
    matrix caaaaug(7, valarray<double>(7));
    matrix guaaaac(7, valarray<double>(7));
    matrix gaaaauc(7, valarray<double>(7));
    matrix cgaaacg(7, valarray<double>(7));
    matrix ggaaacc(7, valarray<double>(7));
    matrix gcaaagc(7, valarray<double>(7));
    matrix aauuuuu(7, valarray<double>(7));
    matrix aauuuuuuu(9, valarray<double>(9));
    matrix uauuuuuua(9, valarray<double>(9));
    matrix ggccccccc(9, valarray<double>(9));
    matrix cgccccccg(9, valarray<double>(9));

    aaaau[0][4] = 0.0037968717806414576;
    gaaac[0][4] = 0.028074937645023688;
    gcccc[0][4] = 0.028074937645023688;

    aacccuu[0][5] = 0.0037525400689852172;
    aacccuu[1][6] = 0.0037525400689852172;
    aacccuu[1][5] = 0.004170771515031873;
    aacccuu[0][6] = 0.004170771515031873;

    aucccau[0][6] = 0.004137135860339378;
    aucccau[1][5] = 0.004137135860339378;

    uacccua[0][6] = 0.0040820768577708028;
    uacccua[1][5] = 0.0040820768577708028;

    cuaaaag[0][6] = 0.029794236432450671;
    cuaaaag[1][5] = 0.00555911755114197;

    caaaaug[0][6] = 0.029511575282733005;
    caaaaug[1][5] = 0.0052693957068543811;

    guaaaac[0][6] = 0.029272684795909918;
    guaaaac[1][5] = 0.0050245378885999878;

    gaaaauc[0][6] = 0.029070803822275754;
    gaaaauc[1][5] = 0.0048176140574817472;

    cgaaacg[0][6] = 0.034144242966177643;
    cgaaacg[1][5] = 0.034144242966177643;

    ggaaacc[0][5] = 0.025748057387594135;
    ggaaacc[1][6] = 0.025748057387594135;
    ggaaacc[0][6] = 0.031384933327228394;
    ggaaacc[1][5] = 0.031384933327228394;

    gcaaagc[0][6] = 0.032195626936620536;
    gcaaagc[1][5] = 0.032195626936620536;

    aauuuuu[0][4] = 0.0037385111560737021;
    aauuuuu[0][5] = 0.0037385111560737021;
    aauuuuu[1][6] = 0.0037385111560737021;
    aauuuuu[0][6] = 0.0041551790391934915;
    aauuuuu[1][5] = 0.0041551790391934915;

    aauuuuuuu[0][4] = 0.00368025229617621;
    aauuuuuuu[0][5] = 0.00368025229617621;
    aauuuuuuu[1][8] = 0.00368025229617621;
    aauuuuuuu[0][6] = 0.0040904270608290343;
    aauuuuuuu[1][7] = 0.0040904270608290343;
    aauuuuuuu[0][7] = 0.0041044537645134214;
    aauuuuuuu[1][6] = 0.0041044537645134214;
    aauuuuuuu[0][8] = 0.0041184804681978085;
    aauuuuuuu[1][5] = 0.0041184804681978085;
    
    uauuuuuua[0][8] = 0.0040345049489403668;
    uauuuuuua[1][5] = 0.0037252446968115789;
    uauuuuuua[1][6] = 0.0037252446968115789;
    uauuuuuua[1][7] = 0.0040062163960097501;
    uauuuuuua[2][8] = 0.0037111004203462705;
    uauuuuuua[3][8] = 0.0037111004203462705;
    uauuuuuua[4][8] = 0.0037111004203462705;

    ggccccccc[0][4] = 0.022541618016111854;
    ggccccccc[0][5] = 0.022541618016111854;
    ggccccccc[1][8] = 0.022541618016111854;
    ggccccccc[0][6] = 0.02747652639862419;
    ggccccccc[1][7] = 0.02747652639862419;
    ggccccccc[0][7] = 0.028127661496104728;
    ggccccccc[1][6] = 0.028127661496104728;
    ggccccccc[0][8] = 0.028778796593585263;
    ggccccccc[1][5] = 0.028778796593585263;

    cgccccccg[0][8] = 0.031403823397438947;
    cgccccccg[1][5] = 0.024534725338365041;
    cgccccccg[1][6] = 0.024534725338365041;
    cgccccccg[1][7] = 0.030026201629414186;
    cgccccccg[2][8] = 0.023845914454352659;
    cgccccccg[3][8] = 0.023845914454352659;
    cgccccccg[4][8] = 0.023845914454352659;

    matrix truth[] = {aaaaa, auaaa, aaaau, gaaac, gcccc, aacccuu, aucccau, uacccua, cuaaaag, caaaaug, guaaaac, gaaaauc,
                        cgaaacg, ggaaacc, gcaaagc, aauuuuu, aauuuuuuu, uauuuuuua, ggccccccc, cgccccccg};

    for (int i = 0; i < 20; i++) {
        cout << seqs[i] << endl;
        RNA thor(seqs[i], false, ener, true, false);
        for (int j = 0; j < seqs[i].length(); j++) {
            cout << j << endl;
            REQUIRE( vect_almost_equal( thor.get_bpp_full()[j], truth[i][j]) );
        }
    }
}
*/
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
