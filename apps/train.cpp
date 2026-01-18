#include <iostream>
#include <cmath>

#include "../hmm/hmm_io.hpp"
#include "../utils/structs_consts_functions.hpp"
#include "../algorithms/baum_welch.hpp"


/**
 * @brief Treniranje skrivenog Markovljevog modela pomoću Baum–Welch algoritma.
 *
 * Ova funkcija provodi iterativno treniranje HMM-a za jedan kromosom:
 * - Učitava prethodno istrenirane parametre ako postoje,
 *   u suprotnom koristi inicijalne parametre.
 * - Učitava trening sekvencu za trenutačni kromosom.
 * - Izvodi više iteracija Baum–Welch algoritma
 *   do konvergencije ili maksimalnog broja iteracija.
 * - Sprema ažurirane HMM parametre na disk.
 *
 * Indeks kromosoma se čuva unutar HMM strukture
 * i koristi se za sekvencijalno treniranje više kromosoma.
 *
 * @note Funkcija je predviđena za uzastopno pozivanje
 *       nad više kromosoma.
 */
int main() {
    HMM hmm;
    
    if (ifstream("../output/trained_hmm_params.txt")) {
        hmm = load_hmm("../output/trained_hmm_params.txt");
    } else {
        hmm = load_hmm("../output/init_hmm_params.txt");
    }

    ifstream in("../output/" + to_string(hmm.chromosome) + "_train_chr.txt");
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


    double prev_ll = -1e100;
    for (int iter = 0; iter < 10; iter++) {
        double ll = 0.0;
        baum_welch_iteration(O, hmm, ll);

        cout << "Iter " << iter << " logL = " << ll << endl;

        if (fabs(ll - prev_ll) < 1e-5) break;
        prev_ll = ll;
    }

    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}