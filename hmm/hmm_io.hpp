#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "../utils/structs_consts_functions.hpp"

using namespace std;


/**
 * @brief Učitava HMM parametre iz datoteke
 * 
 * @param filename Ime datoteke
 * @return HMM Struktura sa učitanim parametrima
 */
HMM load_hmm(const string &filename);


/**
 * @brief Sprema HMM parametre u datoteku
 * 
 * @param hmm HMM parametri
 * @param filename Ime datoteke
 */
void save_hmm(const HMM& hmm, const string& filename);