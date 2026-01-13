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

    double BB, BC, CC, CB;
    compute_transition_probabilities(coords, background.size(), BB, BC, CC, CB);

    hmm.A[0][0] = BB;
    hmm.A[0][1] = BC;
    hmm.A[1][1] = CC;
    hmm.A[1][0] = CB;

    save_hmm(hmm, "../output/init_hmm_params.txt");
    cout << "Inicijalni HMM parametri spremljeni u init_hmm_params.txt\n";
    return 0;
}