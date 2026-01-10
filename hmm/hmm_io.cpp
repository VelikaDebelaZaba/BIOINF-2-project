#include "./hmm_io.hpp"
#include "./hmm.hpp"

using namespace std;


HMM load_hmm(const string &filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Ne mogu otvoriti HMM parametre\n";
        exit(1);
    }

    HMM hmm;
    string tmp;

    if (filename == "../output/trained_hmm_params.txt") {
        in >> tmp >> hmm.chromosome;
    } else {
        hmm.chromosome = -1;
    }

    in >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;   
    in >> hmm.B[0][0] >> hmm.B[0][1] >> hmm.B[0][2] >> hmm.B[0][3];

    in >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
    in >> hmm.B[1][0] >> hmm.B[1][1] >> hmm.B[1][2] >> hmm.B[1][3];

    in >> tmp;
    in >> tmp >> hmm.A[0][0];
    in >> tmp >> hmm.A[0][1];
    in >> tmp >> hmm.A[1][1];
    in >> tmp >> hmm.A[1][0];

    hmm.pi[0] = 0.9;
    hmm.pi[1] = 0.1;

    return hmm;
}


void save_hmm(const HMM& hmm, const string& filename) {
    ofstream out(filename);

    out << "Kromosom: " << hmm.chromosome + 1 << "\n";
    out << "Emisije B (A C G T): ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[0][k] << " ";
    out << "\nEmisije C (A C G T): ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[1][k] << " ";

    out << "\nTranzicije:\n";
    out << "B->B " << hmm.A[0][0] << "\n";
    out << "B->C " << hmm.A[0][1] << "\n";
    out << "C->C " << hmm.A[1][1] << "\n";
    out << "C->B " << hmm.A[1][0] << "\n";
}