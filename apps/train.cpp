#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>

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

// Provjera preklapanja intervala [a,b] s CpG listom (komprimirani indeksi, 1-based)
static bool overlaps_any(int a, int b, const std::vector<CpgRegion>& cpg) {
    // cpg mora biti sortiran po start
    auto it = std::lower_bound(
        cpg.begin(), cpg.end(), a,
        [](const CpgRegion& r, int val){ return r.end < val; }
    );

    if (it != cpg.end()) {
        const auto& r = *it;
        if (!(b < r.start || a > r.end)) return true;
    }
    if (it != cpg.begin()) {
        const auto& r = *(it - 1);
        if (!(b < r.start || a > r.end)) return true;
    }
    return false;
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

    hmm.pi[0] = 0.9;
    hmm.pi[1] = 0.1;

    if (std::ifstream("../output/trained_hmm_params.txt")) {
        hmm = load_hmm("../output/trained_hmm_params.txt");
    } else {
        hmm = load_hmm("../output/init_hmm_params.txt");
    }

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

    // ------- SAMPLING PARAMETRI -------
    const int FLANK = 200;
    const int MIN_POS_WIN = 200;
    const long long MAX_POS_BASES = 5'000'000;
    const int NEG_RATIO = 5;

    // ------- pozitivni prozori -------
    std::vector<std::pair<int,int>> pos_windows;
    long long pos_bases = 0;

    for (const auto& r : coords_chr) {
        int a = std::max(1, r.start - FLANK);
        int b = std::min((int)s.size(), r.end + FLANK);
        if (b - a + 1 < MIN_POS_WIN) continue;

        pos_windows.push_back({a, b});
        pos_bases += (b - a + 1);

        if (pos_bases >= MAX_POS_BASES) break;
    }

    if (pos_windows.empty()) {
        std::cerr << "Nema pozitivnih prozora za chr " << chr << " (mapping coords nije uspio?)" << std::endl;
        return 1;
    }

    // ------- negativni prozori -------
    long long target_neg_bases = (long long)NEG_RATIO * pos_bases;

    std::mt19937 rng(12345 + chr);
    std::uniform_int_distribution<int> len_dist(200, 2000);
    std::uniform_int_distribution<int> start_dist(1, (int)s.size());

    std::vector<std::pair<int,int>> neg_windows;
    long long neg_bases = 0;

    int tries = 0;
    while (neg_bases < target_neg_bases && tries < 2'000'000) {
        tries++;
        int L = len_dist(rng);
        int a = start_dist(rng);
        int b = a + L - 1;
        if (b > (int)s.size()) continue;

        if (overlaps_any(a, b, coords_chr)) continue;

        neg_windows.push_back({a, b});
        neg_bases += L;
    }

    std::cout << "Sampling chr" << chr
              << " pos_bases=" << pos_bases
              << " neg_bases=" << neg_bases
              << " ratio~1:" << (neg_bases / (double)pos_bases) << "\n";

    // ------- složi training sekvence (dinukleotidi) -------
    std::vector<std::vector<int>> sequences;
    sequences.reserve(pos_windows.size() + neg_windows.size());

    for (auto [a,b] : pos_windows) {
        auto O = to_obs_dinuc(s, a, b);
        if ((int)O.size() >= 2) sequences.push_back(std::move(O));
    }
    for (auto [a,b] : neg_windows) {
        auto O = to_obs_dinuc(s, a, b);
        if ((int)O.size() >= 2) sequences.push_back(std::move(O));
    }

    std::shuffle(sequences.begin(), sequences.end(), rng);

    if (sequences.empty()) {
        std::cerr << "Nema validnih training sekvenci (dinukleotidi) za chr " << chr << std::endl;
        return 1;
    }

    // ------- Baum-Welch na mini-sekvencama -------
    double prev_ll = -1e100;
    for (int iter = 0; iter < 10; iter++) {
        double ll = 0.0;
        baum_welch_iteration_multi(sequences, hmm, ll);

        std::cout << "Iter " << iter << " logL = " << ll << std::endl;

        if (fabs(ll - prev_ll) < 1e-3) break;
        prev_ll = ll;
    }

    // sljedeći kromosom (train 1..16)
    if (hmm.chromosome < 16) hmm.chromosome++;
    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
