#include "./viterbij.hpp"
#include <cmath>

static inline double safe_log(double x) {
    if (x <= 0.0 || !std::isfinite(x)) return -1e300;
    return std::log(x);
}

vector<int> viterbi(const vector<int>& O, const HMM& hmm) {
    const int T = (int)O.size();
    if (T <= 0) return {};

    vector<array<double, NSTATE>> delta(T);
    vector<array<int, NSTATE>> prev(T);

    // init t=0
    for (int i = 0; i < NSTATE; i++) {
        delta[0][i] = safe_log(hmm.pi[i]) + safe_log(hmm.B[i][O[0]]);
        prev[0][i] = 0;
    }

    // recursion t=1..T-1
    for (int t = 1; t < T; t++) {
        for (int i = 0; i < NSTATE; i++) {
            double best = -1e300;
            int best_state = 0;

            for (int j = 0; j < NSTATE; j++) {
                const double val = delta[t - 1][j] + safe_log(hmm.A[j][i]);
                if (val > best) {
                    best = val;
                    best_state = j;
                }
            }

            delta[t][i] = best + safe_log(hmm.B[i][O[t]]);
            prev[t][i] = best_state;
        }
    }

    // backtrace
    vector<int> path(T);
    path[T - 1] = (delta[T - 1][0] > delta[T - 1][1]) ? 0 : 1;

    for (int t = T - 2; t >= 0; t--) {
        path[t] = prev[t + 1][path[t + 1]];
    }

    return path;
}
