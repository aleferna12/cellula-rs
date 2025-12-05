//
// Created by aleferna on 31/01/24.
//

#include <iostream>
#include <vector>
#include "sweep_genomes.h"

using namespace sweep;
using namespace std;


void writeGenomeCSV(vector<Genome> &genomes, string &outputfile) {
    ofstream file(outputfile);
    if (not file)
        throw runtime_error("Failed to open file");

    file << GENOME_HEADERS << "\n";
    for (auto &genome : genomes) {
        file << genome.stringRepresentation() << "\n";
    }
}


int main(int argc, char *argv[]) {
    auto genomes = readGenomes(argv[1]);
    string outputfile = argv[2];
    int replicas = stoi(argv[3]);
    int generations = stoi(argv[4]);
    double mut_rate = stod(argv[5]), mustd = stod(argv[6]);
    vector<Genome> mut_genomes;
    for (const auto &genome : genomes) {
        mut_genomes.push_back(genome);
        for (int r = 0; r < replicas; r++) {
            auto mut_genome = genome;
            for (int g = 0; g < generations; g++) {
                mut_genome.MutateGenome(mut_rate, mustd);
                mut_genomes.push_back(mut_genome);  // Copy genome to vector
            }
        }
    }
    writeGenomeCSV(mut_genomes, outputfile);
    return EXIT_SUCCESS;
}
