#include "hmm/hmm_io.hpp"
#include "utils/structs_consts_functions.hpp"
#include "algorithms/viterbij.hpp"
#include "postprocesing/decoded_postprocesing.hpp"
#include "evaluation/evaluation.hpp"



int main() {
    HMM hmm = load_hmm("../output/trained_hmm_params.txt");

    ifstream in("../output/" + to_string(hmm.chromosome) + "_test_chr.txt");
    string s;
    getline(in, s);

    vector<int> O;
    for (char c : s) O.push_back(sym_index(c));
    cout << "Učitana sekvenca za kromosom " << hmm.chromosome << " dužine " << O.size() << endl;

    vector<int> states = viterbi(O, hmm);

    vector<CpgRegion> predicted;
    extract_cpg_islands(predicted, states);
    move_predicted_based_on_lowercase(predicted, hmm.chromosome);
    filter_lenght_and_merge_close_islands(predicted);

    vector<CpgRegion> true_islands;
    load_true_islands(true_islands, hmm.chromosome);
    island_based_evaluation(predicted, true_islands);
    base_pair_evaluation(predicted, true_islands);

    // koristimo funkciju kako bi samo podigli trenutni kromosom za 1 
    save_hmm(hmm, "../output/trained_hmm_params.txt");
}