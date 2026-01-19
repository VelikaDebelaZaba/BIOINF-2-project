#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "../hmm/hmm_io.hpp"
#include "../algorithms/viterbij.hpp"
#include "../postprocesing/decoded_postprocesing.hpp"
#include "../evaluation/evaluation.hpp"
#include "../utils/structs_consts_functions.hpp"

// Pomoć: zadrži samo islands koji se preklapaju sa "keep" intervalom u baznim koordinatama,
// i skrati ih na taj interval (da ne uzimamo rubove prozora)
static void keep_and_clip(std::vector<CpgRegion>& islands, int keep_start_bp, int keep_end_bp) {
    std::vector<CpgRegion> out;
    out.reserve(islands.size());

    for (auto r : islands) {
        int s = std::max(r.start, keep_start_bp);
        int e = std::min(r.end, keep_end_bp);
        if (e >= s) out.push_back({s, e, 0});
    }

    islands.swap(out);
}

int main() {
    HMM hmm = load_hmm("../output/trained_hmm_params.txt");

    hmm.pi[0] = 0.96;
    hmm.pi[1] = 0.04;

    // test kromosomi 17-22
    if (hmm.chromosome < 17) hmm.chromosome = 17;

    std::ifstream in("../output/" + std::to_string(hmm.chromosome) + "_test_chr.txt");
    if (!in) {
        std::cerr << "Ne mogu otvoriti test fajl za kromosom " << hmm.chromosome << "\n";
        return 1;
    }

    std::string s;
    std::getline(in, s);
    if ((int)s.size() < 2) {
        std::cerr << "Premala/prazna sekvenca za kromosom " << hmm.chromosome << "\n";
        return 1;
    }

    // Dinukleotidna opažanja za cijeli kromosom
    std::vector<int> O;
    O.reserve(s.size() - 1);
    for (size_t i = 1; i < s.size(); i++) {
        int x = di_index(s[i - 1], s[i]);
        if (x != -1) O.push_back(x);
    }

    if (O.empty()) {
        std::cerr << "Nema validnih dinukleotida u sekvenci za kromosom " << hmm.chromosome << "\n";
        return 1;
    }

    std::cout << "Učitana sekvenca za kromosom " << hmm.chromosome
              << " (baze=" << s.size()
              << ", dinukleotidi=" << O.size() << ")\n";

    // Windowing parametri (dinukleotidi)
    const int WIN = 5'000'000;       // broj dinukleotida po prozoru
    const int OVERLAP = 50'000;      // preklapanje (dinukleotidi)
    const int STEP = WIN - OVERLAP;  // pomak prozora

    std::vector<CpgRegion> predicted_all;
    predicted_all.reserve(20000);

    int T = (int)O.size(); // dinukleotidi
    for (int start_d = 0; start_d < T; start_d += STEP) {
        int end_d = std::min(start_d + WIN, T);
        int len_d = end_d - start_d;
        if (len_d < 2) break;

        // segment opažanja
        std::vector<int> Oseg;
        Oseg.assign(O.begin() + start_d, O.begin() + end_d);

        // viterbi na segmentu
        std::vector<int> states = viterbi(Oseg, hmm);

        // islands u BAZNIM koordinatama segmenta (zbog extract_cpg_islands mapiranja end+1)
        std::vector<CpgRegion> islands_seg;
        extract_cpg_islands(islands_seg, states);

        // Globalni offset:
        // start_d je 0-based dinukleotid index u globalnom O.
        // Dinukleotid t (0-based) pokriva baze (t+1, t+2) u globalnim bazama.
        // Naš extract vraća bazne koordinate 1-based unutar segmenta.
        // Zato global shift u bazama = start_d
        // (jer segment baza start = start_d+1)
        int base_shift = start_d;

        for (auto& r : islands_seg) {
            r.start += base_shift;
            r.end   += base_shift;
        }

        // Zadrži samo “sredinu” prozora da rubovi ne rade artefakte:
        // keep interval u bazama:
        // - lijevi rub: za prvi prozor zadržimo od početka
        // - inače odrežemo prvih OVERLAP/2 dinukleotida (~ baza)
        // - desni rub: za zadnji prozor zadržimo do kraja
        int keep_left_d  = (start_d == 0) ? start_d : start_d + OVERLAP / 2;
        int keep_right_d = (end_d == T)   ? end_d  : end_d  - OVERLAP / 2;

        // dinukleotidi -> baze: [d_left, d_right] dinukleotidi pokrivaju baze [d_left+1, d_right+1]
        int keep_start_bp = keep_left_d + 1;
        int keep_end_bp   = keep_right_d + 1;

        keep_and_clip(islands_seg, keep_start_bp, keep_end_bp);
        filter_lenght_and_merge_close_islands(islands_seg);

        // dodaj u globalnu listu
        predicted_all.insert(predicted_all.end(), islands_seg.begin(), islands_seg.end());

        std::cerr << "Window " << start_d << "-" << end_d
                  << " (bp keep " << keep_start_bp << "-" << keep_end_bp
                  << "), islands=" << islands_seg.size() << "\n";
    }

    // Sredi/merge
    move_predicted_based_on_lowercase(predicted_all, hmm.chromosome);

    // Truth + eval
    std::vector<CpgRegion> true_islands;
    load_true_islands(true_islands, hmm.chromosome);
    island_based_evaluation(predicted_all, true_islands);
    base_pair_evaluation(predicted_all, true_islands);

    if (hmm.chromosome < 22) hmm.chromosome++;
    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
