#include "./hmm_params_init.hpp"


void hmm_params_init() {
    vector<string> cpg = load_sequences("../../output/clean_positive.txt");
    string background = load_background("../../output/clean_background.txt");
    vector<CpgRegion> coords = load_coords("../../output/coords.txt");

    HMM hmm;

    compute_emission_pos(cpg, hmm.B[1]);      
    compute_emission_bg(background, hmm.B[0]); 

    double BB, BC, CC, CB;
    compute_transition_probabilities(coords, background.size(), BB, BC, CC, CB);

    hmm.A[0][0] = BB;
    hmm.A[0][1] = BC;
    hmm.A[1][1] = CC;
    hmm.A[1][0] = CB;

    ofstream out("../../output/init_hmm_params.txt");
    out << std::fixed << std::setprecision(8);  // format ispisa
    
    out << "Emisije B (A C G T): "
        << hmm.B[0][0] << " "
        << hmm.B[0][1] << " "
        << hmm.B[0][2] << " "
        << hmm.B[0][3] << "\n";

    out << "Emisije C (A C G T): "
        << hmm.B[1][0] << " "
        << hmm.B[1][1] << " "
        << hmm.B[1][2] << " "
        << hmm.B[1][3] << "\n";

    out << "Tranzicije:\n";
    out << "B->B: " << hmm.A[0][0] << "\n";
    out << "B->C: " << hmm.A[0][1] << "\n";
    out << "C->C: " << hmm.A[1][1] << "\n";
    out << "C->B: " << hmm.A[1][0] << "\n";

    cout << "Inicijalni HMM parametri spremljeni u init_hmm_params.txt\n";
    return;
}