#include "../hmm/hmm.hpp"
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
void hmm_params_init();