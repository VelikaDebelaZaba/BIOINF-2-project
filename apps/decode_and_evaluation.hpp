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
void decode_and_evaluate();