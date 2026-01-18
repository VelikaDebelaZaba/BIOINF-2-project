#include "./hmm_io.hpp"
#include "./hmm.hpp"

#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

static bool starts_with(const string& s, const string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

static vector<double> parse_doubles_after_colon(const string& line) {
    size_t pos = line.find(':');
    string tail = (pos == string::npos) ? line : line.substr(pos + 1);

    istringstream iss(tail);
    vector<double> vals;
    double x;
    while (iss >> x) vals.push_back(x);
    return vals;
}

HMM load_hmm(const string &filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Ne mogu otvoriti HMM parametre\n";
        exit(1);
    }

    HMM hmm;

    // default vrijednosti (ako nešto fali u fajlu)
    hmm.chromosome = 1;
    hmm.pi[0] = 0.9;
    hmm.pi[1] = 0.1;

    // sigurni default za A i B (u slučaju da nešto ne učitamo)
    for (int i = 0; i < NSTATE; i++) {
        for (int j = 0; j < NSTATE; j++) hmm.A[i][j] = (i == j) ? 0.99 : 0.01;
        for (int k = 0; k < NSYM; k++) hmm.B[i][k] = 1.0 / NSYM;
    }

    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;

        if (starts_with(line, "Kromosom:")) {
            // format: "Kromosom: <broj>"
            istringstream iss(line);
            string tmp;
            iss >> tmp; // "Kromosom:"
            iss >> hmm.chromosome;
            continue;
        }

        if (starts_with(line, "Emisije B")) {
            auto vals = parse_doubles_after_colon(line);
            if ((int)vals.size() >= NSYM) {
                for (int k = 0; k < NSYM; k++) hmm.B[0][k] = vals[k];
            }
            continue;
        }

        if (starts_with(line, "Emisije C")) {
            auto vals = parse_doubles_after_colon(line);
            if ((int)vals.size() >= NSYM) {
                for (int k = 0; k < NSYM; k++) hmm.B[1][k] = vals[k];
            }
            continue;
        }

        if (starts_with(line, "B->B")) {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> hmm.A[0][0];
            continue;
        }
        if (starts_with(line, "B->C")) {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> hmm.A[0][1];
            continue;
        }
        if (starts_with(line, "C->C")) {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> hmm.A[1][1];
            continue;
        }
        if (starts_with(line, "C->B")) {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> hmm.A[1][0];
            continue;
        }

        if (starts_with(line, "PI:")) {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> hmm.pi[0] >> hmm.pi[1];
            continue;
        }
    }

    return hmm;
}

void save_hmm(const HMM& hmm, const string& filename) {
    ofstream out(filename);
    out << fixed << setprecision(8);

    // zadržavamo isto ponašanje: u trained file spremamo kromosom
    if (filename == "../output/trained_hmm_params.txt") {
        out << "Kromosom: " << hmm.chromosome << "\n";
    }

    // Emisije: ispisujemo NSYM brojeva
    out << "Emisije B: ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[0][k] << " ";
    out << "\n";

    out << "Emisije C: ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[1][k] << " ";
    out << "\n";

    out << "Tranzicije:\n";
    out << "B->B " << hmm.A[0][0] << "\n";
    out << "B->C " << hmm.A[0][1] << "\n";
    out << "C->C " << hmm.A[1][1] << "\n";
    out << "C->B " << hmm.A[1][0] << "\n";

    // PI uvijek spremimo (da bude konzistentno)
    out << "PI: " << hmm.pi[0] << " " << hmm.pi[1] << "\n";
}
