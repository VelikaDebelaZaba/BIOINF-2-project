#pragma once

#include "utils/structs_consts_functions.hpp"
#include "algorithms/forward_backward.hpp"

using namespace std;

/** 
 * @brief Jedna iteracija Baum-Welch algoritma za ažuriranje HMM parametara.
 * 
 * @param O Niz opažanja
 * @param hmm Referenca na HMM strukturu za ažuriranje parametara
 * @param ll Referenca za pohranu log-vjerojatnosti sekvenci
 *
 * @return double Ažurirana log-vjerojatnost sekvenci
*/
double baum_welch_iteration(const vector<int>& O, HMM& hmm, double& ll);
