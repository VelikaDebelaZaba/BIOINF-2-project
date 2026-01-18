#include "./hmm_io.hpp"
#include "./hmm.hpp"
#include <sstream>

using namespace std;


HMM load_hmm(const string &filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Ne mogu otvoriti HMM parametre\n";
        exit(1);
    }

    HMM hmm;
    string line;

    if (filename == "../output/trained_hmm_params.txt") {
        getline(in, line);
        auto pos = line.find(':');
        hmm.chromosome = (pos != string::npos) ? stoi(line.substr(pos + 1)) : 1;
    } else {
        hmm.chromosome = 1;
    }

    auto read_emissions = [&](double emit[NSYM]) {
        string line;
        while (getline(in, line)) {
            if (!line.empty()) break;
        }
        auto pos = line.find(':');
        if (pos == string::npos) {
            cerr << "Neispravan format emisija\n";
            exit(1);
        }
        istringstream iss(line.substr(pos + 1));
        for (int k = 0; k < NSYM; k++) {
            if (!(iss >> emit[k])) {
                cerr << "Neispravan broj emisija\n";
                exit(1);
            }
        }
    };


    read_emissions(hmm.B[0]);
    read_emissions(hmm.B[1]);

    while (getline(in, line)) {
        if (!line.empty()) break;
    }


    if (line.find("Tranzicije") == string::npos) {
        cerr << "Neispravan format tranzicija\n";
        exit(1);
    }

    auto read_transition = [&](double &value) {
        if (!getline(in, line)) {
            cerr << "Neispravan format tranzicija\n";
            exit(1);
        }
        istringstream iss(line);
        string label;
        if (!(iss >> label >> value)) {
            cerr << "Neispravan format tranzicija\n";
            exit(1);
        }
    };

    read_transition(hmm.A[0][0]);
    read_transition(hmm.A[0][1]);
    read_transition(hmm.A[1][1]);
    read_transition(hmm.A[1][0]);

    if (getline(in, line)) {
        auto pos = line.find(':');
        if (pos != string::npos) {
            istringstream iss(line.substr(pos + 1));
            iss >> hmm.pi[0] >> hmm.pi[1];
        }
    }

    if (filename != "../output/trained_hmm_params.txt") {
        hmm.pi[0] = 0.9;
        hmm.pi[1] = 0.1;
    }

    return hmm;
}


void save_hmm(const HMM& hmm, const string& filename) {
    ofstream out(filename);
    out << fixed << setprecision(8);  // format ispisa

    if (filename == "../output/trained_hmm_params.txt") 
        out << "Kromosom: " << hmm.chromosome + 1 << "\n";
    
    out << "Emisije B (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT): ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[0][k] << " ";
    out << "\nEmisije C (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT): ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[1][k] << " ";

    out << "\nTranzicije:\n";
    out << "B->B " << hmm.A[0][0] << "\n";
    out << "B->C " << hmm.A[0][1] << "\n";
    out << "C->C " << hmm.A[1][1] << "\n";
    out << "C->B " << hmm.A[1][0] << "\n";

    if (filename == "../output/trained_hmm_params.txt") {
        out << "PI: " << hmm.pi[0] << " " << hmm.pi[1] << "\n";
    } else {
        out << "PI: 0.9 0.1\n";
    }
}