#ifndef GeneHeader
#define GeneHeader

#include <vector>

class Gene {
public:

    int genetype; //input=0, reg=1, out=2
    int genenr;
    int Boolstate;
    int NewBool; //for updating purposes
    double threshold; //at which value of inputs is this gene "on"

    std::vector<double> w_innode; //input weights coming from input nodes
    std::vector<double> w_regnode; //input weights coming from regulation nodes

    Gene(int type, int nr, int innr, int regnr);

    ~Gene();

    Gene(const Gene &obj); //copy constructor

    //gene functions
    void inline EndCycle() { Boolstate = NewBool; }

    void Mutate(double mu, double mustd);

//private:
    //double mu=0.1;
    //double mustd=0.1;

};


#endif
