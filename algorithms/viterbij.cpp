#include "./viterbij.hpp"


vector<int> viterbi(const vector<int>& O, const HMM& hmm) {
    int T = O.size();

    vector<array<double, NSTATE>> delta(T);
    vector<array<int, NSTATE>> prev(T);

    // inicijalizacija
    for (int i = 0; i < NSTATE; i++) {
        delta[0][i] = log(hmm.pi[i]) + log(hmm.B[i][O[0]]);
        prev[0][i] = 0;
    }

    for (int t = 1; t < T; t++) {
        for (int i = 0; i < NSTATE; i++) {
            double max_val = -1e300;
            int max_state = 0;

            for (int j = 0; j < NSTATE; j++) {
                double val = delta[t-1][j] + log(hmm.A[j][i]);
                if (val > max_val) {
                    max_val = val;
                    max_state = j;
                }
            }

            delta[t][i] = max_val + log(hmm.B[i][O[t]]);
            prev[t][i] = max_state;
        }
    }

    vector<int> path(T);
    path[T-1] = (delta[T-1][0] > delta[T-1][1]) ? 0 : 1;

    for (int t = T-2; t >= 0; t--) {
        path[t] = prev[t+1][path[t+1]];
    }

    return path;
}