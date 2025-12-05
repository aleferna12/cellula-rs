//
// Created by aleferna on 31/01/24.
//

#ifndef REGULATION_EVOLUTION_SWEEP_GENOMES_H
#define REGULATION_EVOLUTION_SWEEP_GENOMES_H

#include <vector>
#include "genome.h"
#include "misc.h"

namespace sweep {

const string GENOME_HEADERS = "innr,regnr,outnr,in_scale_list,reg_threshold_list,reg_w_innode_list,"
                              "reg_w_regnode_list,out_threshold_list,out_w_regnode_list";
const string OUT_HEADERS = "chem,foodp,tau,jkey_dec,jlock_dec,reg_states,out_states,id";
using Input = std::pair<double, double>;
using Output = std::pair<int, std::pair<int, int>>;

struct GenomeTableEntry {
    double chem;
    double foodp;
    int tau;
    int jkey_dec;
    int jlock_dec;
    vector<bool> reg_states;  // We could also convert to a decimal but i think this is easier
    vector<bool> out_states;  // This could be retrieved with jkey and jlock but its more convenient this way
    int id;
};

vector<Genome> readGenomes(const string &inputfile) {
    ifstream file(inputfile);
    if (not file)
        throw runtime_error("Failed to open file");

    string line;
    getline(file, line);
    if (line != GENOME_HEADERS)
        throw runtime_error("Inadequate file headers for genome input. Should be: " + GENOME_HEADERS);

    vector<Genome> genomes{};
    while (getline(file, line)) {
        auto attrs = stringToVector<string>(line, ',');
        auto it = attrs.begin();
        Genome genome;
        Genome::readGenomeInfo(it, genome);
        genomes.push_back(std::move(genome));
    }

    return genomes;
}

void writeGenomeTable(vector<GenomeTableEntry> &genome_table, string &outputfile) {
    ofstream file(outputfile);
    if (not file)
        throw runtime_error("Failed to open file");

    file << OUT_HEADERS << endl;
    for (auto &row : genome_table) {
        file << row.chem << ',' << row.foodp << ',' << row.tau << ',' << row.jkey_dec << ',' << row.jlock_dec << ','
             // We use a ' ' delimiter to prevent pandas from interpreting this field as an integer which is annoying
             << vectorToString(row.reg_states, ' ') << ','
             << vectorToString(row.out_states, ' ') << ','
             << row.id << endl;
    }
}

vector<Input> makeInputs(
        double min_chem,
        double max_chem,
        double step_chem,
        double min_foodparc,
        double max_foodparc,
        double step_foodparc
) {
    vector<Input> inputs {};
    // Don't use <= comparisons for the loops because of double arithmetic imprecision
    for (double chem = min_chem; chem < max_chem; chem += step_chem)
        for (double food = min_foodparc; food < max_foodparc; food += step_foodparc)
            inputs.emplace_back(chem, food);
    return inputs;
}

vector<bool> getGeneStates(vector<Gene> &genes) {
    vector<bool> states;
    states.reserve(genes.size());
    for (auto &gene : genes) {
        states.push_back(gene.Boolstate);
    }
    return states;
}

Output getOutput(Genome &genome, Input &input, int MCSs) {
    for (int i = 0; i < MCSs; ++i) {
        genome.UpdateGeneExpression(array<double, 2>{input.first, input.second}, true);
        genome.FinishUpdate();
    }
    // Tau here is slightly inaccurate because we are not accounting for the fact cells take time to enter dividing state
    return {genome.outputnodes[0].Boolstate + 1, genome.calculateJdecs()};
}

vector<GenomeTableEntry>
makeGenomeTable(vector<Genome> &genomes, vector<Input> &inputs, int MCSs, bool reset) {
    vector<GenomeTableEntry> genome_table {};
    for (size_t i = 0; i < genomes.size(); ++i) {
        auto genome = genomes[i];
        for (auto &input: inputs) {
            if (reset)
                genome.ResetGenomeState();
            auto output = getOutput(genome, input, MCSs);

            GenomeTableEntry row{
                    input.first,
                    input.second,
                    output.first,
                    output.second.first,
                    output.second.second,
                    getGeneStates(genome.regnodes),
                    getGeneStates(genome.outputnodes),
                    (int) i
            };
            genome_table.push_back(row);
        }
    }
    return genome_table;
}

}

#endif //REGULATION_EVOLUTION_SWEEP_GENOMES_H
