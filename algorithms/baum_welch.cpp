#include <cmath>

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

        // zaštita od dijeljenja s nulom
        if (xi_nazivnik == 0.0) xi_nazivnik = 1e-300;
        if (gamma_nazivnik == 0.0) gamma_nazivnik = 1e-300;

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

            // emisije za opažanje u trenutku t
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

    // Emisije za zadnji trenutak (t = T-1) nisu obuhvaćene u petlji gore, pa ih dodamo ovdje.
    {
        int t = T - 1;
        double gamma_den = 0.0;
        for (int i = 0; i < NSTATE; i++)
            gamma_den += alpha[t][i] * beta[t][i];

        if (gamma_den == 0.0) gamma_den = 1e-300;

        for (int i = 0; i < NSTATE; i++) {
            double gamma = (alpha[t][i] * beta[t][i]) / gamma_den;
            B_den[i] += gamma;
            B_num[i][O[t]] += gamma;
        }
    }

/*
    for (int i = 0; i < NSTATE; i++) {
        // zaštita: ako je stanje praktički nedostižno
        if (A_den[i] == 0.0) A_den[i] = 1e-300;
        if (B_den[i] == 0.0) B_den[i] = 1e-300;

        for (int j = 0; j < NSTATE; j++)
            hmm.A[i][j] = A_num[i][j] / A_den[i];

        for (int k = 0; k < NSYM; k++)
            hmm.B[i][k] = B_num[i][k] / B_den[i];
    }
*/
    return ll;
}

double baum_welch_iteration_multi(const vector<vector<int>>& sequences, HMM& hmm, double& ll) {
    double A_num[NSTATE][NSTATE] = {{0}};
    double A_den[NSTATE] = {0};

    double B_num[NSTATE][NSYM] = {{0}};
    double B_den[NSTATE] = {0};

    int used_sequences = 0;

    for (const auto& O : sequences) {
        int T = (int)O.size();
        if (T < 2) continue;

        vector<array<double, NSTATE>> alpha, beta;
        vector<double> c;

        double lseq = forward_scaled(O, hmm, alpha, c);
        if (!std::isfinite(lseq)) continue;     // ako forward pukne, preskoči sekvencu
        ll += lseq;

        backward_scaled(O, hmm, c, beta);

        // --- t = 0..T-2 ---
        for (int t = 0; t < T - 1; t++) {
            double xi_den = 0.0;
            for (int i = 0; i < NSTATE; i++)
                for (int j = 0; j < NSTATE; j++)
                    xi_den += alpha[t][i] * hmm.A[i][j] * hmm.B[j][O[t+1]] * beta[t+1][j];

            double gamma_den = 0.0;
            for (int i = 0; i < NSTATE; i++)
                gamma_den += alpha[t][i] * beta[t][i];

            if (!std::isfinite(xi_den) || xi_den <= 0.0) xi_den = 1e-300;
            if (!std::isfinite(gamma_den) || gamma_den <= 0.0) gamma_den = 1e-300;

            for (int i = 0; i < NSTATE; i++) {
                double gamma = (alpha[t][i] * beta[t][i]) / gamma_den;
                if (!std::isfinite(gamma) || gamma < 0.0) continue;

                A_den[i] += gamma;
                B_den[i] += gamma;
                B_num[i][O[t]] += gamma;

                for (int j = 0; j < NSTATE; j++) {
                    double xi = alpha[t][i] * hmm.A[i][j] * hmm.B[j][O[t+1]] * beta[t+1][j] / xi_den;
                    if (!std::isfinite(xi) || xi < 0.0) continue;
                    A_num[i][j] += xi;
                }
            }
        }

        // --- zadnji simbol za emisije (t = T-1) ---
        {
            int t = T - 1;
            double gamma_den = 0.0;
            for (int i = 0; i < NSTATE; i++) gamma_den += alpha[t][i] * beta[t][i];
            if (!std::isfinite(gamma_den) || gamma_den <= 0.0) gamma_den = 1e-300;

            for (int i = 0; i < NSTATE; i++) {
                double gamma = (alpha[t][i] * beta[t][i]) / gamma_den;
                if (!std::isfinite(gamma) || gamma < 0.0) continue;

                B_den[i] += gamma;
                B_num[i][O[t]] += gamma;
            }
        }

        used_sequences++;
    }

    // --- update pi = prosjek preko sekvenci ---
    // Ako nijedna sekvenca nije korištena, ne diraj parametre (izbjegni NaN)
    if (used_sequences == 0) {
        return ll;
    }

    // --- update A i B s pseudocountom (sprječava nule) ---
    const double A_PSEUDO = 1e-3;
    const double B_PSEUDO = 1e-2;
    const double B_FLOOR = 1e-6;

    for (int i = 0; i < NSTATE; i++) {
        if (!std::isfinite(A_den[i]) || A_den[i] <= 0.0) A_den[i] = 1e-300;
        if (!std::isfinite(B_den[i]) || B_den[i] <= 0.0) B_den[i] = 1e-300;

        double Aden = A_den[i] + A_PSEUDO * NSTATE;
        double Bden = B_den[i] + B_PSEUDO * NSYM;

        /*
        for (int j = 0; j < NSTATE; j++) {
            double num = A_num[i][j];
            if (!std::isfinite(num) || num < 0.0) num = 0.0;
            hmm.A[i][j] = (num + A_PSEUDO) / Aden;
        }
        */
        double bsum = 0.0;
        for (int k = 0; k < NSYM; k++) {
            double num = B_num[i][k];
            if (!std::isfinite(num) || num < 0.0) num = 0.0;
            hmm.B[i][k] = (num + B_PSEUDO) / Bden;
            if (hmm.B[i][k] < B_FLOOR) hmm.B[i][k] = B_FLOOR;
                bsum += hmm.B[i][k];
            }
            if (bsum > 0.0) {
                for (int k = 0; k < NSYM; k++) {
                    hmm.B[i][k] /= bsum;
                }
            }
        }


    return ll;
}
