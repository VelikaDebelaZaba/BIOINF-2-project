#include "../preprocesing/genome_preprocesing.hpp"


int main() {
    vector<CpgRegion> coords;
    string background;
    long chromosome_length = 0;

    const string output_dir = "../output";
    const string genome_path = "../data/ncbi_dataset/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna";
    const int NUM_CHROMOSOMES = 22;

    vector<string> positive_cpg = load_positive_cpg("../data/test.txt", coords);
    background = load_background(genome_path, coords);

    vector<ofstream> chromosome_out_files;
    ofstream out1, out2, coords_out;
    open_output_files(NUM_CHROMOSOMES, output_dir, chromosome_out_files, out1, out2, coords_out);


    for (int chr = 1; chr <= NUM_CHROMOSOMES; chr++) {
        vector<lowerCaseRegions> lowercaseCoords;
        string chromosome = load_chromosome(genome_path, lowercaseCoords, chr);
        cout << "Duzina kromosoma " << chr << ": " << chromosome.size() << endl;
        chromosome_length += chromosome.size();
        chromosome_out_files[chr - 1] << chromosome << "\n";
        for (const auto &l : lowercaseCoords) 
            chromosome_out_files[chr - 1] << l.start << " " << l.end << "\n";
    }

    cout << "Broj pozitivnih CpG otoka: " << positive_cpg.size() << endl;
    cout << "Duzina originalnog kromosoma: " << chromosome_length << endl;
    cout << "Duzina backgrounda nakon ciscenja: " << background.size() << endl;

    for (const auto &s : positive_cpg) out1 << s << "\n";
    out2 << background;
    for (const auto &c : coords) coords_out << c.chromosome << " " << c.start << " " << c.end << "\n";

    for (int chr = 0; chr < NUM_CHROMOSOMES; chr++) chromosome_out_files[chr].close();
    out1.close();
    out2.close();
    coords_out.close();

    return 0;
}