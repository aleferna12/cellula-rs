//
// Created by aleferna on 23-02-2023.
//
#include <iostream>
#include <vector>
#include <algorithm>
#include "sweep_genomes.h"

using namespace std;
using namespace sweep;

int main(int argc, char *argv[]) {
    vector<string> args(argv + 1, argv + argc);
    for (auto &arg : args) {
        if (arg == "-h" or arg == "--help") {
            cout << "Usage: sweep_genome <inputfile> <outputfile> <min_chem> "
                    "<max_chem> <step_chem> <min_foodparc> <max_foodparc> <step_foodparc> [MCSs] [reset]"
                    "\n\n"
                    "Where:\n"
                    "\t-'inputfile' must be a CSV file similar to the ones used to backup "
                    "cell data but only containing genome attributes "
                    "(" << GENOME_HEADERS << ")\n"
                    "\t-'MCSs' controls for how many time-steps the simulated genome will receive the same inputs "
                    "(default: 50, type: INT)\n"
                    "\t-'reset' determines whether the genome is reset before a new combination of inputs "
                    "is given (default: true, type: BOOL)\n"
                    "\t-Food parameters should be given in terms of division parcels: "
                    "food parcels = food / (grn_update_period * (divtime + divdur) / metabperiod)\n"
                    "\t-All numeric arguments besides 'MCSs' can be doubles or integers\n"
                    "\t-The input range is inclusive in the interval [min_val, max_val]"
                    "\n\n"
                    "Description: Creates a CSV file with the output of a genome for all "
                    "combinations of inputs given by the input parameters of the program" << endl;
            return EXIT_SUCCESS;
        }
    }
    if (args.size() < 8) {
        cerr << "Inadequate arguments, try: sweep_genome -h" << endl;
        return EXIT_FAILURE;
    }
    int MCSs = 50;
    if (args.size() > 8) {
        MCSs = stoi(args[8]);
        if (MCSs < 1) {
            cerr << "'MCSs' argument must be 1 or higher" << endl;
            return EXIT_FAILURE;
        }
    }
    bool reset = true;
    if (args.size() > 9) {
        reset = args[9] != "false" or args[9] != "0";
    }

    auto genomes = readGenomes(args[0]);
    auto inputs = makeInputs(
        stod(args[2]),
        stod(args[3]),
        stod(args[4]),
        stod(args[5]),
        stod(args[6]),
        stod(args[7])
    );
    auto genome_table = makeGenomeTable(genomes, inputs, MCSs, reset);
    writeGenomeTable(genome_table, args[1]);
    return EXIT_SUCCESS;
}