#include "./forward_backward.hpp"    


double forward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
) {
    int T = O.size();
    alpha.assign(T, {});
    c.assign(T, 0.0);

    // t = 0
    for (int i = 0; i < NSTATE; i++) {
        alpha[0][i] = hmm.pi[i] * hmm.B[i][O[0]];
        c[0] += alpha[0][i];
    }

    c[0] = 1.0 / c[0];
    for (int i = 0; i < NSTATE; i++)
        alpha[0][i] *= c[0];

    // t = 1..T-1
    for (int t = 1; t < T; t++) {
        for (int j = 0; j < NSTATE; j++) {
            alpha[t][j] = 0.0;
            for (int i = 0; i < NSTATE; i++)
                alpha[t][j] += alpha[t-1][i] * hmm.A[i][j];
            alpha[t][j] *= hmm.B[j][O[t]];
            c[t] += alpha[t][j];
        }

        c[t] = 1.0 / c[t];
        for (int j = 0; j < NSTATE; j++)
            alpha[t][j] *= c[t];
    }

    double loglik = 0.0;
    for (int t = 0; t < T; t++)
        loglik -= log(c[t]);

    return loglik;
}


void backward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
) {
    int T = O.size();
    beta.assign(T, {});

    for (int i = 0; i < NSTATE; i++)
        beta[T-1][i] = c[T-1];

    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < NSTATE; i++) {
            beta[t][i] = 0.0;
            for (int j = 0; j < NSTATE; j++)
                beta[t][i] += hmm.A[i][j] * hmm.B[j][O[t+1]] * beta[t+1][j];
            beta[t][i] *= c[t];
        }
    }
}