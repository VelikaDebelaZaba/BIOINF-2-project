#include "./baum_welch.hpp"
//#include <iostream>


double baum_welch_iteration(const vector<int>& O, HMM& hmm, double& ll) {
    double A_num[NSTATE][NSTATE] = {{0}};
    double A_den[NSTATE] = {0};

    double B_num[NSTATE][NSYM] = {{0}};
    double B_den[NSTATE] = {0};

    int T = O.size();

    vector<array<double, NSTATE>> alpha, beta;
    vector<double> c;

    ll += forward_scaled(O, hmm, alpha, c);
    backward_scaled(O, hmm, c, beta);

    for (int t = 0; t < T-1; t++) {

        double xi_nazivnik = 0.0;
        for (int i = 0; i < NSTATE; i++)
            for (int j = 0; j < NSTATE; j++)
                xi_nazivnik +=
                    alpha[t][i] *
                    hmm.A[i][j] *
                    hmm.B[j][O[t+1]] *
                    beta[t+1][j];

        double gamma_nazivnik = 0.0;
        for (int i = 0; i < NSTATE; i++)
            gamma_nazivnik += alpha[t][i] * beta[t][i];

        for (int i = 0; i < NSTATE; i++) {
            double gamma = (alpha[t][i] * beta[t][i]) / gamma_nazivnik;

            // azuriranje pi vrijednosti u HMM modelu
            // PROVJERITI treba li se skalirati jer pri ispisu dosta pada
            // kaže sin kada smanjimo veličinu uzorka da bi trebalo biti uredu
            /*if (t == 0) {
                cout << gamma << "\n";
                hmm.pi[i] = gamma;
            }*/

            A_den[i] += gamma;
            B_den[i] += gamma;
            B_num[i][O[t]] += gamma;

            for (int j = 0; j < NSTATE; j++) {
                double xi =
                    alpha[t][i] *
                    hmm.A[i][j] *
                    hmm.B[j][O[t+1]] *
                    beta[t+1][j] / xi_nazivnik;

                A_num[i][j] += xi;
            }
        }
    }

    for (int i = 0; i < NSTATE; i++) {
        for (int j = 0; j < NSTATE; j++)
            hmm.A[i][j] = A_num[i][j] / A_den[i];

        for (int k = 0; k < NSYM; k++)
            hmm.B[i][k] = B_num[i][k] / B_den[i];
    }

    return ll;
}