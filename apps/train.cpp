#include "./train.hpp"


void train_hmm() {
    HMM hmm;
    
    if (ifstream("../../output/trained_hmm_params.txt")) {
        hmm = load_hmm("../../output/trained_hmm_params.txt");
    } else {
        hmm = load_hmm("../../output/init_hmm_params.txt");
        hmm.chromosome = 1;
    }

    ifstream in("../../output/" + to_string(hmm.chromosome) + "_train_chr.txt");
    string s;
    getline(in, s);

    vector<int> O;
    for (char c : s) O.push_back(sym_index(c));
    cout << "Učitana sekvenca za kromosom " << hmm.chromosome << " dužine " << O.size() << endl;


    double prev_ll = -1e100;
    for (int iter = 0; iter < 10; iter++) {
        double ll = 0.0;
        baum_welch_iteration(O, hmm, ll);

        cout << "Iter " << iter << " logL = " << ll << endl;

        if (fabs(ll - prev_ll) < 1e-5) break;
        prev_ll = ll;
    }

    save_hmm(hmm, "../../output/trained_hmm_params.txt");
}