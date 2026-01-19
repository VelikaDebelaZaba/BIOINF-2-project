#include <iostream>
#include <cmath>
#include <algorithm>


#include "../hmm/hmm_io.hpp"
#include "../hmm/hmm.hpp"
#include "../utils/structs_consts_functions.hpp"
#include "../algorithms/baum_welch.hpp"

// Ucitaj sekvencu + lowercase intervale iz train_chr fajla
static void load_chr_and_lowercase(const std::string& filename, std::string& seq, std::vector<lowerCaseRegions>& lc) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Ne mogu otvoriti: " << filename << std::endl;
        exit(1);
    }

    std::getline(in, seq);

    lc.clear();
    int s, e;
    while (in >> s >> e) {
        lc.push_back({s, e});
    }
    std::sort(lc.begin(), lc.end(), [](auto& a, auto& b){ return a.start < b.start; });
}

// Broj lowercase baza striktno prije pozicije x (1-based) u originalnom koordinatnom sustavu
static int lowercase_before(const std::vector<lowerCaseRegions>& lc, int x) {
    long long removed = 0;
    for (const auto& r : lc) {
        if (r.start >= x) break;
        int end = std::min(r.end, x - 1);
        if (end >= r.start) removed += (end - r.start + 1);
    }
    return (int)removed;
}

// Map original 1-based coord -> komprimirana 1-based coord (uppercase-only)
static int map_orig_to_comp(const std::vector<lowerCaseRegions>& lc, int orig_pos) {
    int removed = lowercase_before(lc, orig_pos);
    return orig_pos - removed;
}

// Dinukleotidna opažanja za prozor [start1, end1] (1-based indeksi u s)
// O će imati duljinu (end1-start1) jer su to parovi baza.
static std::vector<int> to_obs_dinuc(const std::string& s, int start1, int end1) {
    std::vector<int> O;

    int L = end1 - start1 + 1;
    if (L < 2) return O;

    O.reserve(L - 1);

    // parovi: (start1, start1+1) ... (end1-1, end1)
    for (int pos = start1 + 1; pos <= end1; pos++) {
        char prev = s[(size_t)(pos - 2)]; // baza na (pos-1) u 1-based
        char cur  = s[(size_t)(pos - 1)]; // baza na pos u 1-based
        int x = di_index(prev, cur);
        if (x != -1) O.push_back(x);
    }

    return O;
}

int main() {
    HMM hmm;

    if (std::ifstream("../output/trained_hmm_params.txt")) {
        hmm = load_hmm("../output/trained_hmm_params.txt");
    } else {
        hmm = load_hmm("../output/init_hmm_params.txt");
    }

    hmm.pi[0] = 0.9;
    hmm.pi[1] = 0.1;

    const int chr = hmm.chromosome;

    // ucitaj sekvencu + lowercase intervale
    std::string s;
    std::vector<lowerCaseRegions> lc;
    load_chr_and_lowercase("../output/" + std::to_string(chr) + "_train_chr.txt", s, lc);

    std::cout << "Učitana sekvenca za kromosom " << chr << " dužine " << s.size() << std::endl;

    // ucitaj coords i filtriraj na ovaj chr (originalne koordinate)
    std::vector<CpgRegion> all = load_coords("../output/coords.txt");
    std::vector<CpgRegion> coords_chr_orig;
    coords_chr_orig.reserve(all.size());
    for (const auto& r : all) if (r.chromosome == chr) coords_chr_orig.push_back(r);

    // mapiraj coords u komprimirani indeks (uppercase-only)
    std::vector<CpgRegion> coords_chr;
    coords_chr.reserve(coords_chr_orig.size());
    for (const auto& r : coords_chr_orig) {
        int cs = map_orig_to_comp(lc, r.start);
        int ce = map_orig_to_comp(lc, r.end);
        if (cs < 1) cs = 1;
        if (ce > (int)s.size()) ce = (int)s.size();
        if (cs <= ce) coords_chr.push_back({cs, ce, chr});
    }
    std::sort(coords_chr.begin(), coords_chr.end(), [](auto& a, auto& b){ return a.start < b.start; });

    // ------- SEMI-SUPERVIZIJA: maska dozvoljenih stanja -------
    // CpG regije su "clamped" na stanje 1, a udaljene regije na stanje 0.
    const int NEG_MARGIN = 200;
    std::vector<char> base_is_cpg(s.size() + 1, 0); // 1-based
    std::vector<char> base_near_cpg(s.size() + 1, 0);

    for (const auto& r : coords_chr) {
        int a = std::max(1, r.start);
        int b = std::min((int)s.size(), r.end);
        for (int pos = a; pos <= b; pos++) base_is_cpg[pos] = 1;
        int na = std::max(1, r.start - NEG_MARGIN);
        int nb = std::min((int)s.size(), r.end + NEG_MARGIN);
        for (int pos = na; pos <= nb; pos++) base_near_cpg[pos] = 1;
    }

    const int T_full = (int)s.size() - 1; // broj dinukleotida
    const int CHUNK_D = 1'000'000; // dinukleotidi po chunku
   
    std::vector<std::vector<int>> sequences;
    std::vector<std::vector<std::array<double, NSTATE>>> masks;

    for (int start_d = 0; start_d < T_full; start_d += CHUNK_D) {
        int end_d = std::min(start_d + CHUNK_D, T_full);
        int start_bp = start_d + 1;
        int end_bp = end_d + 1;

        auto O = to_obs_dinuc(s, start_bp, end_bp);
        if ((int)O.size() < 2) continue;

        std::vector<std::array<double, NSTATE>> mask;
        mask.reserve(O.size());

        for (int d = start_d; d < end_d; d++) {
            int b1 = d + 1;
            int b2 = d + 2;
            if (b2 > (int)s.size()) break;

            if (base_is_cpg[b1] || base_is_cpg[b2]) {
                mask.push_back({0.0, 1.0});
            } else if (!base_near_cpg[b1] && !base_near_cpg[b2]) {
                mask.push_back({1.0, 0.0});
            } else {
                mask.push_back({1.0, 1.0});
            }
        }

        if (mask.size() == O.size()) {
            sequences.push_back(std::move(O));
            masks.push_back(std::move(mask));
        }
    }

    if (sequences.empty()) {
        std::cerr << "Nema validnih training sekvenci (dinukleotidi) za chr " << chr << std::endl;
        return 1;
    }

    // ------- Baum-Welch na mini-sekvencama -------
    double prev_ll = -1e100;
    for (int iter = 0; iter < 10; iter++) {
        double ll = 0.0;
        baum_welch_iteration_multi_masked(sequences, masks, hmm, ll);

        std::cout << "Iter " << iter << " logL = " << ll << std::endl;

        if (fabs(ll - prev_ll) < 1e-3) break;
        prev_ll = ll;
    }

    // sljedeći kromosom (train 1..16)
    if (hmm.chromosome < 16) hmm.chromosome++;
    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
