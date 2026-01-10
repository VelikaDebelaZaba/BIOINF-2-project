#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "../utils/structs_consts_functions.hpp"

using namespace std;

/**
 * @brief Viterbi algoritam za pronalaženje najvjerojatnije sekvence stanja za dani kromosom.
 * 
 * @param O Sekvenca opažanja
 * @param hmm HMM parametri
 * 
 * @return vector<int> Najvjerojatnija sekvenca stanja za dani kromosom
 */
vector<int> viterbi(const vector<int>& O, const HMM& hmm);
