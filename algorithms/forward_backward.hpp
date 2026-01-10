#pragma once

#include <vector>
#include <array>
#include <cmath>

#include "utils/structs_consts_functions.hpp"

using namespace std;


/**
 * @brief Forward algoritam sa skaliranjem kako bi izbjegli padanje
 * vrijednosti alpha u 0.
 * 
 * @param O Sekvenca opažanja
 * @param hmm HMM parametri
 * @param alpha Matrica za pohranu forward varijabli
 * @param c Vektor skalirajućih faktora
 * 
 * @return double Log-vjerojatnost sekvence
 */
double forward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
);

/**
 * @brief Backward algoritam sa skaliranjem kako bi izbjegli padanje
 * vrijednosti beta u 0. Koristi skalirajuće faktore iz forward algoritma.
 * 
 * @param O Niz opažanja
 * @param hmm HMM parametri
 * @param c Vektor skalirajućih faktora iz forward algoritma
 * @param beta Matrica za pohranu backward varijabli
 */
void backward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
);
