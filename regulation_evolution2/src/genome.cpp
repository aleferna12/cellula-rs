#include "genome.h"
#include "random.h"
#include "output.h"
#include "misc.h"

//default constructor required for init in other class
Genome::Genome() {
    innr = 0;
    regnr = 0;
    outnr = 0;
}

Genome::Genome(int in, int reg, int out, vector<double> const &inscale) {
    InitGenome(in, reg, out, inscale);
}

void Genome::InitGenome(int in, int reg, int out, vector<double> const &inscale) {
    int i;

    innr = in;
    regnr = reg;
    outnr = out;

    // Inputs must match in order of magnitude as to not introduce bias in the effect of mutations
    inputscale = inscale;

    //for (i=0; i<in; i++)
    //{
    //  inputnodes.push_back(*(new Gene(0, i, innr, regnr)));
    //}

    for (i = 0; i < reg; i++) {
        //regnodes.push_back(*(new Gene(1, i+innr, innr, regnr)));
        regnodes.emplace_back(1, i + innr, innr, regnr);
    }

    for (i = 0; i < out; i++) {
        //outputnodes.push_back(*(new Gene(2, i+innr+regnr, innr, regnr)));
        outputnodes.emplace_back(2, i + innr + regnr, innr, regnr);
    }
}

Genome::~Genome() {
    //cout<<"destructed"<<endl;
}

//copy constructor
Genome::Genome(const Genome &Parent) {
    innr = Parent.innr;
    regnr = Parent.regnr;
    outnr = Parent.outnr;

    inputscale = Parent.inputscale;
    //inputnodes=Parent.inputnodes;
    regnodes = Parent.regnodes;
    outputnodes = Parent.outputnodes;
    for (auto &n: regnodes) {
        n.Boolstate = 0;
        n.NewBool = 0;
    }
    for (auto &n: outputnodes) {
        n.Boolstate = 0;
        n.NewBool = 0;
    }
}

void Genome::ResetGenomeState() {
    for (auto &n: regnodes) {
        n.Boolstate = 0;
        n.NewBool = 0;
    }
    for (auto &n: outputnodes) {
        n.Boolstate = 0;
        n.NewBool = 0;
    }
}

void Genome::ReadFromFile(char *filename) {
    std::ifstream ifs;
    string line;
    int i, nodetype, nodenr;
    double scale, thresh;

    ifs.open(filename, std::ifstream::in);

    if (ifs.is_open()) {

        getline(ifs, line);
        stringstream strstr(line);

        //first read the nr of nodes
        strstr >> innr >> regnr >> outnr;

        for (i = 0; i < regnr; i++) {
            regnodes.push_back(*(new Gene(1, i + innr, innr, regnr)));
        }

        for (i = 0; i < outnr; i++) {
            outputnodes.push_back(*(new Gene(2, i + innr + regnr, innr, regnr)));
        }

        getline(ifs, line);
        strstr.clear();
        strstr.str(std::string());
        strstr << line;
        for (i = 0; i < innr; i++) {
            strstr >> scale;
            inputscale.push_back(scale);
        }

        //now read all the interaction strengths
        getline(ifs, line);
        while (line.length()) {
            strstr.clear();
            strstr.str(std::string());
            strstr << line;
            //stringstream strstr(line);
            //read the straightforward cell variables from the line
            strstr >> nodetype >> nodenr >> thresh;

            //if node is regulatory
            if (nodetype == 1) {
                regnodes[nodenr - innr].threshold = thresh;
                for (auto &w: regnodes[nodenr - innr].w_innode) {
                    strstr >> w;
                }
                for (auto &w: regnodes[nodenr - innr].w_regnode) {
                    strstr >> w;
                }
            }
                //if node is output
            else {
                outputnodes[nodenr - innr - regnr].threshold = thresh;
                for (auto &w: outputnodes[nodenr - innr - regnr].w_regnode) {
                    strstr >> w;
                }
            }
            //get next line
            getline(ifs, line);
        }

    } else {
        cerr << "Genome::ReadFromFile error: could not open file. exiting ..." << endl;
        exit(1);
    }
}

void Genome::WriteToFile(char *filename) {
    std::ofstream ofs;
    ofs.open(filename, std::ofstream::out);

    //first read the nr of nodes
    ofs << innr << " " << regnr << " " << outnr << endl;
    for (auto i: inputscale) {
        ofs << i << " ";
    }
    ofs << endl;
    //now write all the interaction strengths

    for (const auto &n: regnodes) {
        ofs << "1 " << n.genenr << " " << n.threshold;
        for (const auto w: n.w_innode) {
            ofs << " " << w;
        }
        for (const auto w: n.w_regnode) {
            ofs << " " << w;
        }
        ofs << endl;
    }
    for (const auto &n: outputnodes) {
        ofs << "2 " << n.genenr << " " << n.threshold;
        for (const auto w: n.w_regnode) {
            ofs << " " << w;
        }
        ofs << endl;
    }
    ofs.flush();
    ofs.close();
}

void Genome::UpdateGeneExpression(const array<double, 2> &input, bool sync_cells) {
    int i, j;
    double newval;
    vector<double> v_input;

    //transpose input data for normalisation
    for (i = 0; i < input.size(); i++) {
        v_input.push_back(inputscale[i] * (double) input[i]);
    }
    //cerr <<"updating..."<<endl;
    //update the regulatory genes
    for (i = 0; i < regnodes.size(); i++) {
        newval = 0.;
        //regulation by input nodes
        for (j = 0; j < regnodes[i].w_innode.size(); j++) {
            newval += v_input[j] * regnodes[i].w_innode[j];
        }
        //regulation by regulation nodes
        for (j = 0; j < regnodes[i].w_regnode.size(); j++) {
            newval += (double) (regnodes[j].Boolstate) * regnodes[i].w_regnode[j];
        }
        //regulation node is on or off depending on regulation
        if (newval > regnodes[i].threshold) {
            regnodes[i].NewBool = 1;
        } else {
            regnodes[i].NewBool = 0;
        }
    }

    //update output nodes
    for (i = 0; i < outputnodes.size(); i++) {
        newval = 0.;
        for (j = 0; j < outputnodes[i].w_regnode.size(); j++) {
            newval += (double) (regnodes[j].Boolstate) * outputnodes[i].w_regnode[j];
        }

        //output node is on or off depending on regulation
        if (newval > outputnodes[i].threshold) {
            outputnodes[i].NewBool = 1;
        } else {
            outputnodes[i].NewBool = 0;
        }

    }

    //only update BoolState if we update cells asynchronously
    if (!sync_cells) FinishUpdate();
}

//this function puts the NewBool state of nodes into BoolState
void Genome::FinishUpdate() {
    int i;
    //update the actual expression
    for (i = 0; i < regnodes.size(); i++) {
        regnodes[i].EndCycle();
    }
    for (i = 0; i < outputnodes.size(); i++) {
        outputnodes[i].EndCycle();
    }
}

void Genome::MutateGenome(double mu, double mustd) {
    for (auto &n: inputscale) {
        if (RANDOM() < mu) {
            n += RANDNORMAL(0., mustd);
        }
    }

    for (auto & regnode : regnodes) {
        regnode.Mutate(mu, mustd);
        regnode.Boolstate = 0;
        regnode.NewBool = 0;
    }

    for (auto & outputnode : outputnodes) {
        outputnode.Mutate(mu, mustd);
        outputnode.Boolstate = 0;
        outputnode.NewBool = 0;
    }

}

string Genome::stringRepresentation() {
    auto genome_stream = stringstream();
    genome_stream << innr << ',' << regnr << ',' << outnr << ',';
    genome_stream << vectorToString(inputscale, ' ') << ',';

    vector<double> reg_thres {};
    vector<double> reg_w_in {};
    vector<double> reg_w_reg {};
    for (auto &g: regnodes) {
        reg_thres.push_back(g.threshold);
        // Unpack the weight vectors
        reg_w_in.insert(reg_w_in.end(), g.w_innode.begin(), g.w_innode.end());
        reg_w_reg.insert(reg_w_reg.end(), g.w_regnode.begin(), g.w_regnode.end());
    }
    genome_stream << vectorToString(reg_thres, ' ') << ',';
    genome_stream << vectorToString(reg_w_in, ' ') << ',';
    genome_stream << vectorToString(reg_w_reg, ' ') << ',';

    vector<double> out_thres {};
    vector<double> out_w_reg {};
    for (auto &g: outputnodes) {
        out_thres.push_back(g.threshold);
        // Unpack the weight vectors
        out_w_reg.insert(out_w_reg.end(), g.w_regnode.begin(), g.w_regnode.end());
    }
    genome_stream << vectorToString(out_thres, ' ') << ',';
    genome_stream << vectorToString(out_w_reg, ' ');
    return genome_stream.str();
}

vector<string>::iterator Genome::readGenomeInfo(vector<string>::iterator it, Genome &genome) {
    int innr = stoi(*it); ++it;
    int regnr = stoi(*it); ++it;
    int outnr = stoi(*it); ++it;
    genome.innr = innr;
    genome.regnr = regnr;
    genome.outnr = outnr;
    genome.inputscale = stringToVector<double>(*it, ' '); ++it;

    vector<double> reg_thres = stringToVector<double>(*it, ' '); ++it;
    vector<double> reg_w_in = stringToVector<double>(*it, ' '); ++it;
    vector<double> reg_w_reg = stringToVector<double>(*it, ' '); ++it;
    for (int i = 0; i < regnr; ++i) {
        Gene gene = Gene(1, i + innr, innr, regnr);
        gene.threshold = reg_thres[i];
        for (int j = 0; j < innr; ++j)
            gene.w_innode[j] = reg_w_in[i * innr + j];
        for (int j = 0; j < regnr; ++j)
            gene.w_regnode[j] = reg_w_reg[i * regnr + j];
        genome.regnodes.push_back(gene);
    }

    vector<double> out_thres = stringToVector<double>(*it, ' '); ++it;
    vector<double> out_w_reg = stringToVector<double>(*it, ' '); ++it;
    for (int i = 0; i < outnr; ++i) {
        Gene gene = Gene(2, i + innr + regnr, innr, regnr);
        gene.threshold = out_thres[i];
        for (int j = 0; j < regnr; ++j)
            gene.w_regnode[j] = out_w_reg[i * regnr + j];
        genome.outputnodes.push_back(gene);
    }
    return it;
}

pair<int, int> Genome::calculateJdecs() {
    int jkey_dec = 0;
    int jlock_dec = 0;
    int exp_2 = 1;
    int key_lock_len = (outnr - 1) / 2;
    for (int i = 1; i < key_lock_len + 1; ++i) {
        if (outputnodes[i].Boolstate == 1)
            jkey_dec += exp_2;
        if (outputnodes[key_lock_len + i].Boolstate == 1)
            jlock_dec += exp_2;
        exp_2 *= 2;
    }
    return {jkey_dec, jlock_dec};
}

