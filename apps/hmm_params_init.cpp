#include <algorithm>

#include "../hmm/hmm.hpp"
#include "../hmm/hmm_io.hpp"
#include "../utils/structs_consts_functions.hpp"

/**
 * @brief Inicijalizacija parametara skrivenog Markovljevog modela (HMM).
 *
 * Funkcija izračunava početne parametre HMM-a:
 * - Emisijske vjerojatnosti za pozadinsko i CpG stanje
 *   na temelju frekvencija nukleotida.
 * - Prijelazne vjerojatnosti između stanja
 *   na temelju poznatih CpG koordinata i duljine genoma.
 *
 * Inicijalni parametri se spremaju u datoteku `init_hmm_params.txt`
 * i koriste se kao početna točka za Baum–Welch treniranje.
 *
 * @note Pretpostavlja se da je predobrada već izvršena
 *       i da potrebne ulazne datoteke postoje.
 */
int main() {
    vector<string> cpg = load_sequences("../output/clean_positive.txt");
    string background = load_background("../output/clean_background.txt");
    vector<CpgRegion> coords = load_coords("../output/coords.txt");

    HMM hmm;

    compute_emission_pos(cpg, hmm.B[1]);      
    compute_emission_bg(background, hmm.B[0]); 

    // Uzmemo tranzicije iz jednog kromosoma (npr. chr1) i pravu duljinu chr1 iz train fajla
    const int CHR_FOR_INIT = 1;

    // filtriraj coords na odabrani kromosom
    vector<CpgRegion> coords_chr;
    coords_chr.reserve(coords.size());
    for (const auto& r : coords) {
        if (r.chromosome == CHR_FOR_INIT) coords_chr.push_back(r);
    }

    // sortiraj po startu (compute_transition_probabilities to pretpostavlja)
    sort(coords_chr.begin(), coords_chr.end(), [](const CpgRegion& a, const CpgRegion& b) {
        return a.start < b.start;
    });

    // duljina kromosoma = duljina prve linije u "1_train_chr.txt"
    ifstream chr_in("../output/" + to_string(CHR_FOR_INIT) + "_train_chr.txt");
    string chr_seq;
    getline(chr_in, chr_seq);

    double BB, BC, CC, CB;
    compute_transition_probabilities(coords_chr, (int)chr_seq.size(), BB, BC, CC, CB);


    hmm.A[0][0] = BB;
    hmm.A[0][1] = BC;
    hmm.A[1][1] = CC;
    hmm.A[1][0] = CB;

    save_hmm(hmm, "../output/init_hmm_params.txt");
    cout << "Inicijalni HMM parametri spremljeni u init_hmm_params.txt\n";
    return 0;
}