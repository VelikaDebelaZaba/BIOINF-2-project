#include "./hmm.hpp"


vector<string> load_sequences(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti fajl: " << filename << endl;
        exit(1);
    }

    vector<string> seqs;
    string line;
    while (getline(file, line)) {
        if (!line.empty())
            seqs.push_back(line);
    }
    return seqs;
}


string load_background(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti background: " << filename << endl;
        exit(1);
    }

    string genome, line;
    while (getline(file, line)) {
        genome += line;
    }
    return genome;
}


vector<CpgRegion> load_coords(const string &filename) {
    vector<CpgRegion> coords;
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti coords fajl: " << filename << endl;
        exit(1);
    }

    int chr, s, e;
    while (file >> chr >> s >> e) {
        coords.push_back({s, e});
    }
    return coords;
}


void compute_emission_pos(const vector<string> &seqs, double emit[16]) {
    long long counts[NSYM] = {0};
    long long total = 0;

    for (const auto &s : seqs) {
        if (s.size() < 2) continue;
        for (size_t i = 0; i + 1 < s.size(); i++) {
            int idx = dinuc_index(s[i], s[i + 1]);
            if (idx >= 0) {
                counts[idx]++;
                total++;

            }
        }
    }

   for (int k = 0; k < NSYM; k++)
        emit[k] = (total > 0) ? counts[k] / static_cast<double>(total) : 0.0;
}


void compute_emission_bg(const string &bg, double emit[16]) {
    long long counts[NSYM] = {0};
    long long total = 0;


    if (bg.size() >= 2) {
        for (size_t i = 0; i + 1 < bg.size(); i++) {
            int idx = dinuc_index(bg[i], bg[i + 1]);
            if (idx >= 0) {
                counts[idx]++;
                total++;
            }
        }
    }


    for (int k = 0; k < NSYM; k++)
        emit[k] = (total > 0) ? counts[k] / static_cast<double>(total) : 0.0;
}


void compute_transition_probabilities(
    const vector<CpgRegion> &coords,
    int chromosome_length,
    double &p_BB, double &p_BC,
    double &p_CC, double &p_CB
) {
    // --- 1. Prosječna dužina CpG otoka ---   
    double sum_C = 0;
    for (auto &r : coords)
        sum_C += (r.end - r.start + 1);

    double L_C = sum_C / coords.size();

    // --- 2. Prosječna dužina background segmenta ---
    double sum_B = 0;
    int num_B_segments = 0;

    // prije prvog CpG otoka
    if (coords[0].start > 1) {
        sum_B += coords[0].start - 1;
        num_B_segments++;
    }

    // između CpG otoka
    for (size_t i = 1; i < coords.size(); i++) {
        int bg_len = coords[i].start - coords[i - 1].end - 1;
        if (bg_len > 0) {
            sum_B += bg_len;
            num_B_segments++;
        }
    }

    // nakon zadnjeg CpG otoka
    if (coords.back().end < chromosome_length) {
        sum_B += chromosome_length - coords.back().end;
        num_B_segments++;
    }

    if (num_B_segments == 0) {
        cerr << "Greška: nema background segmenata!" << endl;
        exit(1);
    }

    double L_B = sum_B / num_B_segments;


    p_CB = 1.0 / L_C;
    p_CC = 1.0 - p_CB;

    p_BC = 1.0 / L_B;
    p_BB = 1.0 - p_BC;

    cout << "Prosjecna duzina CpG otoka L_C: " << L_C << endl;
    cout << "Prosjecna duzina background segmenta L_B: " << L_B << endl;
}