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
        
        void sum_left_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT);
        void sum_right_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT);
        void sum_exterior_basepairs(int ii, int jj, double q_bound_ij, double exp_neg_gstack_over_RT);
        void gradient_sum_exterior_basepairs(int ii, int jj, double q_bound_ij, double exp_neg_gstack_over_RT);
        void gradient_sum_right_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT);
        void gradient_sum_left_interior_loops(int ii, int jj, double qbij_over_full, double exp_neg_gstack_over_RT);
        void calc_bppS_entry(int ii, int jj);
        void calc_bppGradientS_entry(int ii, int jj);
        void calc_bpp();
        void calc_bpp_gradient();

        bool isCircular;
        bool calcBPP;
        bool calcGrad;
        int nn; // number of bases
        std::string sequence;
        vect energies;
        matrix g_base_pair;
        //partition matrices
        matrix partitionBound;
        matrix partitionS;
        matrix partition;
        tensor gradientBound;
        tensor gradientS;
        tensor gradient;
        //bpp matrices
        matrix bppS;
        matrix bpp;
        tensor bppGradientS;
        tensor bppGradient;
};
