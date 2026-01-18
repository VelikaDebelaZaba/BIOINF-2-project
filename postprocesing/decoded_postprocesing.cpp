#include "./decoded_postprocesing.hpp"


void extract_cpg_islands(vector<CpgRegion>& islands, vector<int>& states) {
    int start_d = -1; // start u dinukleotid indeksima (1-based)

    for (int i = 0; i < (int)states.size(); i++) {
        if (states[i] == 1 && start_d == -1) start_d = i + 1;

        if ((states[i] == 0 || i == (int)states.size() - 1) && start_d != -1) {
            int end_d = (states[i] == 0) ? i : i + 1; // dinukleotid end (1-based)

            // MAPIRANJE: dinukleotid [start_d, end_d] pokriva baze [start_d, end_d+1]
            int start_bp = start_d;
            int end_bp = end_d + 1;

            islands.push_back({start_bp, end_bp, 0}); // chr nije bitan ovdje
            start_d = -1;
        }
    }
}



void move_predicted_based_on_lowercase(vector<CpgRegion>& predicted, int chr_number) {
    vector<CpgRegion> lowercaseCoords;
    ifstream in("../output/" + to_string(chr_number) + "_test_chr.txt");
    string s;
    getline(in, s);

    while(getline(in, s)) {
        size_t space = s.find(' ', 0);
        int start = stoll(s.substr(0, space));
        int end   = stoll(s.substr(space + 1, s.size() - space - 1));
        lowercaseCoords.push_back({start, end});
    }


    for (auto& p : predicted) {
        int offset = 0;

        for (const auto& lc : lowercaseCoords) {
            if (lc.start <= p.start + offset) {
                offset += (lc.end - lc.start + 1);
            } else {
                break;
            }
        }

        p.start += offset;
        p.end += offset;

    }
}


void filter_lenght_and_merge_close_islands(vector<CpgRegion>& islands) {
    if (islands.empty()) return;

    sort(islands.begin(), islands.end(), 
        [](const CpgRegion& a, const CpgRegion& b) {
            return a.start < b.start;
        });

    vector<CpgRegion> temp;
    temp.push_back(islands[0]);

    for (int i = 1; i < islands.size(); i++) {
        if (islands[i].start - temp.back().end <= MERGE_DISTANCE) {
            temp.back().end = max(temp.back().end, islands[i].end);
        } else {
            temp.push_back(islands[i]);
        }
    }

    islands.clear();

    for (const auto& region : temp) {
        if (region.end - region.start + 1 >= MIN_CPG_LEN) {
            islands.push_back(region);
        }
    }
}


void load_true_islands(vector<CpgRegion>& islands, int chr_number) {
    ifstream in("../output/coords.txt");
    int chr, start, end;

    while (in >> chr >> start >> end) {
        if (chr == chr_number) {
            islands.push_back({start, end});
        }
    } 
}



