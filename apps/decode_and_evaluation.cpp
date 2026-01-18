#include "../hmm/hmm_io.hpp"
#include "../algorithms/viterbij.hpp"
#include "../postprocesing/decoded_postprocesing.hpp"
#include "../evaluation/evaluation.hpp"


/**
 * @brief Dekodiranje sekvence i evaluacija predviđenih CpG otoka.
 *
 * Funkcija provodi fazu dekodiranja i evaluacije:
 * - Učitava istrenirane HMM parametre.
 * - Primjenjuje Viterbijev algoritam za određivanje
 *   najvjerojatnijeg niza stanja.
 * - Iz dekodiranih stanja izdvaja predviđene CpG otoke.
 * - Prilagođava predikcije prema lowercase regijama genoma.
 * - Filtrira i spaja bliske ili kratke CpG otoke.
 * - Uspoređuje predviđene otoke s poznatim CpG otocima koristeći:
 *     - evaluaciju na razini otoka
 *     - evaluaciju na razini parova baza
 *
 * Nakon evaluacije, indeks kromosoma se povećava
 * i sprema za sljedeći ciklus dekodiranja.
 *
 * @note Funkcija je namijenjena evaluaciji na testnim kromosomima.
 */ 
int main() {
    HMM hmm = load_hmm("../output/trained_hmm_params.txt");

    ifstream in("../output/" + to_string(hmm.chromosome) + "_test_chr.txt");
    string s;
    getline(in, s);

    vector<int> O;
    if (s.size() >= 2) {
        for (size_t i = 0; i + 1 < s.size(); i++) {
            int idx = dinuc_index(s[i], s[i + 1]);
            if (idx >= 0) O.push_back(idx);
        }
    }

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

    return 0;
}