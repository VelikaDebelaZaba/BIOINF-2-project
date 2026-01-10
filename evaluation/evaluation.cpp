#include "evaluation.hpp"


void island_based_evaluation(const vector<CpgRegion>& predicted, const vector<CpgRegion>& truth) {
    int TP = 0;

    for (const auto& p : predicted) {
        for (const auto& t : truth) {
            int s = max(p.start, t.start);
            int e = min(p.end, t.end);

            if (e >= s) {
                TP++;
                break;
            }
        }
    }

    int FP = predicted.size() - TP;
    int FN = truth.size() - TP;

    double precision = TP / double(TP + FP);
    double recall = TP / double(TP + FN);

    cout << "=== Island-based evaluation ===\n";
    cout << "True Positive = " << TP << "\n";
    cout << "False Positive = " << FP << "\n";
    cout << "False Negative = " << FN << "\n";
    cout << "Precision = " << precision << "\n";
    cout << "Recall = " << recall << "\n";
}


void base_pair_evaluation(const vector<CpgRegion>& predicted, const vector<CpgRegion>& truth) {
    long long pred_len = 0;
    long long truth_len = 0;

    for (const auto& p : predicted)
        pred_len += (p.end - p.start + 1);

    for (const auto& t : truth)
        truth_len += (t.end - t.start + 1);

    long long overlap_len = 0;

    size_t i = 0, j = 0;
    while (i < predicted.size() && j < truth.size()) {
        int s = max(predicted[i].start, truth[j].start);
        int e = min(predicted[i].end,   truth[j].end);

        if (s <= e)
            overlap_len += (e - s + 1);

        if (predicted[i].end < truth[j].end)
            i++;
        else
            j++;
    }

    long long TP = overlap_len;
    long long FP = pred_len  - overlap_len;
    long long FN = truth_len - overlap_len;

    double precision = TP / double(TP + FP);
    double recall    = TP / double(TP + FN);

    cout << "=== Base-pair evaluation ===\n";
    cout << "True Positive = " << TP << "\n";
    cout << "False Positive = " << FP << "\n";
    cout << "False Negative = " << FN << "\n";
    cout << "Precision = " << precision << "\n";
    cout << "Recall = " << recall << "\n";
}