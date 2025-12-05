#include "gene.h"
#include "random.h"

Gene::Gene(int type, int nr, int innr, int regnr) {
    genetype = type;
    genenr = nr;
    //this is a regulatory node
    if (type == 1) {
        //w_innode.resize(innr,0.0); Could init all to 0, but is a bit weird
        for (int i = 0; i < innr; i++) {
            w_innode.push_back(RANDOM() * 2. - 1.); //right version
            //w_innode.push_back(RANDOM()); //test version
        }
    }

    //this is a regulatory or output node
    if (type >= 1) {
        //w_regnode.resize(regnr,0.0);
        for (int i = 0; i < regnr; i++) {
            w_regnode.push_back(RANDOM() * 2. - 1.);//right version
            //w_regnode.push_back(RANDOM()); //test version
        }
    }

    threshold = RANDOM() * 2. - 1.;
    Boolstate = 0;
    NewBool = 0;
}

Gene::Gene(const Gene &obj) {
    genetype = obj.genetype;
    genenr = obj.genenr;
    threshold = obj.threshold;

    w_innode = obj.w_innode;
    w_regnode = obj.w_regnode;

    Boolstate = 0;
    NewBool = 0;

}

Gene::~Gene() {
    //destructor doesn't need to do anything except be there
}


void Gene::Mutate(double mu, double mustd) {
    int i;
    //mutate regulation coming from input nodes
    for (i = 0; i < w_innode.size(); i++) {
        if (RANDOM() < mu) {
            w_innode[i] += RANDNORMAL(0., mustd);
        }
    }

    //mutate regulation coming from regulatory nodes
    for (i = 0; i < w_regnode.size(); i++) {
        if (RANDOM() < mu) {
            w_regnode[i] += RANDNORMAL(0., mustd);
        }
    }

    //mutate activation threshold
    if (RANDOM() < mu) {
        threshold += RANDNORMAL(0., mustd);
    }

}
