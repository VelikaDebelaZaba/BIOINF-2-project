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
void train_hmm();