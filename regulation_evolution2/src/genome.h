#ifndef GenomeHeader
#define GenomeHeader

#include <cstdio>
#include <array>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include "gene.h"

using namespace std;

class Genome {
public:


    int innr, regnr, outnr;

    //transformation of inputs for normalisation
    vector<double> inputscale;

    //Storage of the the GRN state
    //vector < Gene > inputnodes;
    vector<Gene> regnodes;
    vector<Gene> outputnodes;

    //constructor and destructor
    Genome();

    Genome(int in, int reg, int out, vector<double> const &inscale);

    void InitGenome(int in, int reg, int out, const vector<double> &inscale); //constructor body for handy init in classes

    ~Genome();

    Genome(const Genome &Parent); //copy constructor

    //read and write genome
    void ResetGenomeState();

    static vector<string>::iterator readGenomeInfo(vector<string>::iterator it, Genome &genome);

    void ReadFromFile(char *filename);

    void WriteToFile(char *filename);

    //genomic functions
    void UpdateGeneExpression(const array<double, 2> &input, bool sync_cells);

    void FinishUpdate();

    void MutateGenome(double mu, double mustd);

    pair<int, int> calculateJdecs();

    string stringRepresentation();
};


#endif
