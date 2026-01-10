#include "../preprocesing/genome_preprocesing.hpp"


/**
 * @brief Predobrada genomskih podataka i priprema ulaznih datoteka za HMM.
 *
 * Ova funkcija provodi cjelokupni postupak predobrade:
 * - Učitava poznate koordinate pozitivnih CpG otoka.
 * - Učitava i čisti pozadinsku (background) genomsku sekvencu.
 * - Dijeli genom na pojedinačne kromosome.
 * - Bilježi lowercase (maskirane) regije za svaki kromosom.
 * - Sprema per-kromosomske sekvence i pripadne metapodatke.
 *
 * U izlazni direktorij se zapisuju:
 * - Očišćene pozitivne CpG sekvence.
 * - Očišćena pozadinska sekvenca.
 * - Datoteke s kromosomskim sekvencama i lowercase anotacijama.
 * - Datoteka s koordinatama CpG otoka.
 *
 * Funkcija se izvršava jednom prije inicijalizacije i treniranja HMM-a.
 *
 * @note Putanje do datoteka su trenutno zadane u kodu (hardcoded).
 */
void preprocess();