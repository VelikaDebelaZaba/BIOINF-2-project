#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main() {

    const int NUM_TRAIN_CHR = 16;   // 1–16 train
    const int NUM_TOTAL_CHR = 22;   // ukupno kromosoma

    cout << "=== Detekcija CpG otoka pipeline ===\n";
    cout << "\n[1/4] Pretprocesiranje\n";
    system("./bin/preprocess");

    cout << "\n[2/4] Inicijalizacija HMM parametara\n";
    system("./bin/hmm_init");


    cout << "\n[3/4] Treniranje HMM (Baum-Welch)\n";

    for (int chr = 1; chr <= NUM_TRAIN_CHR; chr++) {
        cout << "\n--- Treniranje na kromosomu " << chr << " ---\n";
        system("./bin/train");
    }


    cout << "\n[4/4] Dekodiranje i evaluacija\n";

    for (int chr = NUM_TRAIN_CHR + 1; chr <= NUM_TOTAL_CHR; chr++) {
        cout << "\n--- Dekodiranje i evaluacija na kromosomu " << chr << " ---\n";
        system("./bin/decode_and_evaluation");
    }

    cout << "\n=== Pipeline finished successfully ===\n";
    return 0;
}
