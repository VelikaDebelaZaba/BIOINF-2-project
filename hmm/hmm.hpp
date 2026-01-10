#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "utils/structs_consts_functions.hpp"

using namespace std;


/**
 * @brief Učitava CpG sekvence iz .txt doteteke, svaka linija je jedna sekvenca
 * 
 * @param filename Ime doteteke
 * @return vector<string> Vektor sekvenci
 * 
 * Napomena: Pretpostavlja se da su sekvence već očišćene i spremne za upotrebu.
 */
vector<string> load_sequences(const string &filename);


/**
 * @brief Učitava genom oćišćen od CpG otoka i malih slova iz .txt doteteke
 * 
 * @param filename Ime doteteke
 * @return string Sekvenca pozadinskog genoma
 * 
 * Napomena: Pretpostavlja se da je sekvenca već očišćena i spremna za upotrebu.
 */
string load_background(const string &filename);


/**
 * @brief Učitava koordinate CpG otoka iz doteteke, preskače broj kromosoma.
 * 
 * @param filename Ime doteteke
 * @return vector<CpgRegion> Vektor koordinata CpG otoka
 * 
 * Napomena: Pretpostavlja se da je doteteka ispravno formatirana.
 */
vector<CpgRegion> load_coords(const string &filename);


/**
 * @brief Računa emisijske vjerojatnosti za CpG stanje na temelju pozitivnih sekvenci
 * 
 * @param seqs Vektor CpG sekvenci
 * @param emit Niz od 4 elementa za pohranu vjerojatnosti A,C,G,T
 */
void compute_emission_pos(const vector<string> &seqs, double emit[4]);


/**
 * @brief Računa emisijske vjerojatnosti za pozadinsko stanje na temelju pozadinskog genoma
 * 
 * @param bg Sekvenca pozadinskog genoma
 * @param emit Niz od 4 elementa za pohranu vjerojatnosti A,C,G,T
 */
void compute_emission_bg(const string &bg, double emit[4]);


/**
 * @brief Računa prijelazne vjerojatnosti između stanja na temelju koordinata CpG otoka
 * i ukupne duljine kromosoma
 * 
 * @param coords Vektor koordinata CpG otoka
 * @param chromosome_length Ukupna duljina kromosoma
 * @param p_BB Referenca za pohranu vjerojatnosti B->B
 * @param p_BC Referenca za pohranu vjerojatnosti B->C
 * @param p_CC Referenca za pohranu vjerojatnosti C->C
 * @param p_CB Referenca za pohranu vjerojatnosti C->B
 */
void compute_transition_probabilities(
    const vector<CpgRegion> &coords,
    int chromosome_length,
    double &p_BB, double &p_BC,
    double &p_CC, double &p_CB
);