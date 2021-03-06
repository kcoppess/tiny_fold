// incomplete version using long arrays rather than matrix and tensor
#include <string>
#include <cmath>
#include <valarray>
#include <vector>

typedef std::vector< std::vector< std::valarray<double> > > tensor;
typedef std::vector< std::valarray<double> > matrix;
typedef std::valarray<double> vect;

class RNA {
    public:
        RNA(std::string seq, bool type, vect ener, bool wantBPP, bool grad);
        ~RNA();
        bool is_circular();
        int get_length();
        std::string get_sequence();
        vect get_energy();
        double get_partition();
        vect get_gradient();
        void update_energy(vect ener);
        matrix get_bpp_full();
        double get_bpp(int i, int j);
        tensor get_bpp_gradient_full();
        vect get_bpp_gradient(int i, int j);
    
    private:
        void calc_partition();
        void calc_gBasePair();
        int find_energy_index(int i, int j);
        double hairpin(double gHP);
        double interior(double gBP, char loop);
        void calc_gradient();
        void calc_bpp();
        void calc_bpp_gradient();

        bool isCircular;
        bool calcBPP;
        int nn; // number of bases
        std::string sequence;
        vect energies;
        double* g_base_pair = NULL;
        double* everything = NULL; //partitionBound, partitionS, partition, bppS, bpp
        double* everything_gradient = NULL; //gradientBound, gradientS, gradient, bppGradientS, bppGradient
};
