#include <string>
#include <cmath>
#include <valarray>
#include <vector>

typedef std::vector< std::vector< std::valarray<double> > > tensor;
typedef std::vector< std::valarray<double> > matrix;
typedef std::valarray<double> vect;

class RNA {
    public:
        RNA(std::string seq, bool type, vect ener);
        ~RNA();
        bool is_circular();
        int get_length();
        std::string get_sequence();
        double get_energy();
        double get_partition();
        vect get_gradient();
        void update_energy(vect ener);
    
    private:
        void calc_partition();
        void calc_gBasePair(matrix& gBP);
        double hairpin(double gHP);
        double interior(double gBP, char loop);
        void calc_gradient();

        bool isCircular;
        int nn; // number of bases
        std::string sequence;
        vect energies;
        matrix partitionBound;
        matrix partitionS;
        matrix partition;
        tensor gradientBound;
        tensor gradientS;
        tensor gradient;
};
