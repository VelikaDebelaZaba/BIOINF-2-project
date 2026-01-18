#pragma once


/**
 * Globalne konstante HMM-a
 * 2 stanja: 0 = Background (B), 1 = CpG (C)
 * 16 simbola: dinukleotidi AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
 */
constexpr int NSTATE = 2;
constexpr int NSYM   = 16;


/**
 * Struktura za pohranu HMM parametara
 *
 * A[i][j]  - prijelazna vjerojatnost iz stanja i u stanje j
 * B[i][k]  - emisijska vjerojatnost stanja i za simbol k
 * pi[i]    - inicijalna vjerojatnost stanja i
 * chromosome - broj kromosoma na kojem je HMM zadnje treniran
 */
struct HMM {
    double A[NSTATE][NSTATE];
    double B[NSTATE][NSYM];
    double pi[NSTATE];
    int chromosome;
};


/**
 * Pretvara nukleotid u indeks:
 * A -> 0, C -> 1, G -> 2, T -> 3
 */
inline int sym_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}


/**
 * Pretvara dinukleotid u indeks:
 * AA -> 0, AC -> 1, AG -> 2, AT -> 3,
 * CA -> 4, CC -> 5, CG -> 6, CT -> 7,
 * GA -> 8, GC -> 9, GG -> 10, GT -> 11,
 * TA -> 12, TC -> 13, TG -> 14, TT -> 15
 */
inline int dinuc_index(char first, char second) {
    int first_idx = sym_index(first);
    int second_idx = sym_index(second);

    if (first_idx < 0 || second_idx < 0) return -1;

    return first_idx * 4 + second_idx;
}


/**
 * CpgRegion struktura za pohranu koordinata CpG otoka
 * start      - početna pozicija (1-based)
 * end        - završna pozicija (1-based)
 * chromosome - kroj kromosoma
 */
struct CpgRegion { 
    int start;
    int end;
    int chromosome;
};


/**
 * Struktura za pohranu koordinata malih slova u genomu
 * start - početna pozicija (1-based)
 * end   - završna pozicija (1-based)
 * 
 * Napomena: stavljeno u 1-based radi konzistentnosti sa CpgRegion strukturama
 * koje su fiksno postavljene u 1-based sustav zbog USC koordinata.
 */
struct lowerCaseRegions {
    int start;
    int end;
};
