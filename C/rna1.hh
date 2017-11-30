#include <string>
#include <cmath>
#include <vector>

typedef std::vector< std::vector<double> > matrix;

class RNA {
    public:
        RNA(std::string seq, bool type, std::vector<double> ener);
        ~RNA();
        bool is_circular();
        int get_length();
        std::string get_sequence();
        double get_energy();
        double get_partition();
        std::vector<double> get_gradient();
        void update_energy(std::vector<double> ener);
    
    private:
        void calc_partition();
        void calc_gBasePair(matrix& gBP);
        double hairpin(double gHP);
        double interior(double gBP, char loop);
        void calc_gradient();

        bool isCircular;
        int nn; // number of bases
        std::string sequence;
        std::vector<double> energies;
        matrix partitionBound;
        matrix partition;
        std::vector<matrix> gradientBound;
        std::vector<matrix> gradient;
};
