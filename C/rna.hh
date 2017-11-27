#include <string>
#include <cmath>
#include <vector>
using namespace std;

typedef vector<vector<double>> matrix;

class RNA {
    public:
        RNA(string seq, bool type, vector<double> ener);
        ~RNA();
        bool is_circular();
        int get_length();
        string get_sequence();
        double get_energy();
        double get_partition();
        void update_energy(vector<double> ener);
    private:
        void calc_partitionBound();
        void calc_partition();
        void calc_gBasePair(matrix& gBP);

        double g_loop;
        bool isCircular;
        int nn; // number of bases
        string sequence;
        vector<double> energies;
        matrix partitionBound;
        matrix partition;
};
