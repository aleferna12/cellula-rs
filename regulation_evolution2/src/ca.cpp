/*

Copyright 1995-2006 Roeland Merks, Nick Savill

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/


/* CA.cpp: implementation of Glazier & Graner's Cellular Potts Model */

// This code derives from a Cellular Potts implementation written around 1995
// by Nick Savill

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <bitset>
#include "misc.h"
#include "sticky.h"
#include "random.h"
#include "ca.h"
#include "parameter.h"
#include "dish.h"
#include "crash.h"
#include "hull.h"

#define ZYGFILE(Z) <Z.xpm>
#define XPM(Z) Z ## _xpm
#define ZYGXPM(Z) XPM(Z)

//Leonie's sign function
template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* define default zygote */
/* NOTE: ZYGOTE is normally defined in Makefile!!!!!! */
// #ifndef ZYGOTE
// #define ZYGOTE init
// //#include "init.xpm"
// #else
// #include ZYGFILE(ZYGOTE)
// #endif

/* STATIC DATA MEMBER INITIALISATION */
double copyprob[BOLTZMANN];

const int CellularPotts::nx[25] = {0, 0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 1, 2, 2, 1, -1, -2, -2, -1, 0, 2, 0, -2};
const int CellularPotts::ny[25] = {0, -1, 0, 1, 0, -1, 1, 1, -1, -2, 0, 2, 0, -2, -1, 1, 2, 2, 1, -1, -2, -2, 0, 2, 0};

const int CellularPotts::nbh_level[5] = {0, 4, 8, 20, 24}; //0:self; 1:van Neumann; 2:Moore; etc...

extern Parameter par;


/** PRIVATE **/

using namespace std;

void CellularPotts::BaseInitialisation(vector<Cell> *cells) {
    CopyProb(par.T);
    cell = cells;
    if (par.neighbours >= 1 && par.neighbours <= 4)
        n_nb = nbh_level[par.neighbours];
    else
        throw runtime_error("Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).");

}

CellularPotts::CellularPotts(vector<Cell> *cells,
                             const int sx, const int sy) {

    sigma = nullptr;
    frozen = false;
    thetime = 0;
    zygote_area = 0;

    auto rule = stringToVector<float>(par.key_lock_weights, ' ');
    max_key_lock_dec = int(pow(2, rule.size()));
    KL_strengths = allocateKLStrengths(rule, max_key_lock_dec, par.dynamic_adh);

    BaseInitialisation(cells);
    sizex = sx;
    sizey = sy;

    AllocateSigma(sx, sy);


    // fill borders with special border state
    for (int x = 0; x < sizex; x++) {
        sigma[x][0] = -1;
        sigma[x][sizey - 1] = -1;
    }
    for (int y = 0; y < sizey; y++) {
        sigma[0][y] = -1;
        sigma[sizex - 1][y] = -1;
    }

    if (par.neighbours >= 1 && par.neighbours <= 4)
        n_nb = nbh_level[par.neighbours];
    else
        throw runtime_error("Panic in CellularPotts: parameter neighbours invalid (choose [1-4])");
}

CellularPotts::CellularPotts() {

    sigma = nullptr;
    sizex = 0;
    sizey = 0;
    frozen = false;
    thetime = 0;
    zygote_area = 0;
    KL_strengths = nullptr;
    max_key_lock_dec = 0;

    CopyProb(par.T);

    // fill borders with special border state
    for (int x = 0; x < sizex; x++) {
        sigma[x][0] = -1;
        sigma[x][sizey - 1] = -1;
    }
    for (int y = 0; y < sizey; y++) {
        sigma[0][y] = -1;
        sigma[sizex - 1][y] = -1;
    }
    if (par.neighbours >= 1 && par.neighbours <= 4)
        n_nb = nbh_level[par.neighbours];
    else
        throw runtime_error("Panic in CellularPotts: parameter neighbours invalid (choose [1-4])");
}

// destructor (virtual)
CellularPotts::~CellularPotts() {
    if (sigma) {
        free(sigma[0]);
        free(sigma);
        sigma = nullptr;
    }
    if (KL_strengths) {
        delete[] KL_strengths;
        KL_strengths = nullptr;
    }
}

float * CellularPotts::allocateKLStrengths(vector<float> &rule, int max_kl_dec, bool dynamic_adh) {
    auto *KL_array = new float[max_kl_dec * max_kl_dec] {};
    for (int i = 0; i < max_kl_dec; ++i)
        for (int j = 0; j < max_kl_dec; ++j) {
            int index = i * max_kl_dec + j;
            if (dynamic_adh)
                KL_array[index] = calculateKLStrength(i, j, rule);
        }
    return KL_array;
}

float CellularPotts::calculateKLStrength(unsigned long key_dec, unsigned long  lock_dec, const vector<float> &rule) {
    auto key = bitset<32>{key_dec};
    auto lock = bitset<32>{lock_dec};
    auto matches = ~(key ^ lock);  // XNOR operation
    float res = 0;
    for (int i = 0; i < int(rule.size()); ++i)
        // bitset[0] is the last digit, so we need to iterate in reverse
        if (matches[rule.size() - 1 - i] == 1)
            res += rule[i];
    return res;
}


void CellularPotts::AllocateSigma(int sx, int sy) {

    sizex = sx;
    sizey = sy;

    sigma = (int **) malloc(sizex * sizeof(int *));
    if (sigma == nullptr)
        MemoryWarning();

    sigma[0] = (int *) malloc(sizex * sizey * sizeof(int));
    if (sigma[0] == nullptr)
        MemoryWarning();


    {
        for (int i = 1; i < sizex; i++)
            sigma[i] = sigma[i - 1] + sizey;
    }

    /* Clear CA plane */
    {
        for (int i = 0; i < sizex * sizey; i++)
            sigma[0][i] = 0;
    }

}

//this function is used in ReadBackup in Dish to fill the ca plane
//AND to set the cell's moments (assumes cells have been initialised!)
int CellularPotts::SetNextSigma(int sig) {
    //the plane has a 1 px boundary on all size, therefore we place the pixels
    //within that
    static int xcount = 1, ycount = 1;

    if (xcount >= sizex - 1 || ycount >= sizey - 1) {
        return 1;
    }

    sigma[xcount][ycount] = sig;
    if (sig) {
        (*cell)[sig].area++;
        (*cell)[sig].AddSiteToMoments(xcount, ycount);
    }

    ycount++;
    if (ycount == sizey - 1) {
        ycount = 1;
        xcount++;
    }
    return 0;

}

void CellularPotts::InitializeEdgeList(bool init_single_cell_center) {
    // unordered_set<long long int> edgeSet;
    int neighbour;
    int x, y, xp, yp;
    int c, cp;
    int xstart = 0;
    int xend = sizex - 1;
    int ystart = 0;
    int yend = sizey - 1;

//  bool init_single_cell_center=true;

    //cout<<"edgevector is size "<<edgeSetVector.size_map()<<endl;;

    if (init_single_cell_center && sizex > 200 && sizey > 200) {
        xstart = (sizex - 1) / 2 - 100;
        xend = (sizex - 1) / 2 + 100;
        ystart = (sizey - 1) / 2 - 100;
        yend = (sizey - 1) / 2 + 100;
    }


    for (x = xstart; x < xend; x++) {
        for (y = ystart; y < yend; y++) {
            for (int nb = 1; nb <= n_nb; nb++) {
                c = sigma[x][y]; //c = getSigma here
                xp = nx[nb] + x;  //check neighbour
                yp = ny[nb] + y;

                if (par.periodic_boundaries) {

                    // since we are asynchronic, we cannot just copy the borders once
                    // every MCS

                    if (xp <= 0)
                        xp = sizex - 2 + xp;
                    if (yp <= 0)
                        yp = sizey - 2 + yp;
                    if (xp >= sizex - 1)
                        xp = xp - sizex + 2;
                    if (yp >= sizey - 1)
                        yp = yp - sizey + 2;

                    cp = sigma[xp][yp];
                } else if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
                    cp = -1;
                else
                    cp = sigma[xp][yp]; //if boundaries deal appropriately, else set the neighbour
                //if c != cp and cp is not boundary and sgn(x) is the sign function
                // (check for boundary state in the neigh.)
                if (cp != c && cp != -1 && c != -1 && sgn(cp) + sgn(c) >= 0) {
                    //edgeSetpair.insert({x,y,xp,yp});
                    edgeSetVector.insert({x, y, xp, yp});
                    // cout << edgeSetVector.subset_size() <<endl;
                    // edgeSetVector.push_back({x,y,xp,yp});
                    // edgeSetVector.insert_in_subset(edgeSetVector.subset_size());
                    // cout << edgeSetVector.subset_size() << endl;
                    //edgeSetpair.insert({xp,yp,x,y});
                    edgeSetVector.insert({xp, yp, x, y});
                    // edgeSetVector.push_back({xp,yp,x,y});
                    // edgeSetVector.insert_in_subset(edgeSetVector.subset_size());
                    // edgeVector.push_back({x,y,xp,yp});
                }
            }
        }
    }
}

int CellularPotts::DeltaHWithMedium(int x, int y, PDE *PDEfield) {
    double DH = 0.;
    int i, sxy, sxyp;
    int neighsite;

    /* Compute energydifference *IF* the copying were to occur */
    sxy = sigma[x][y];
    sxyp = MEDIUM; // getSigma[xp][yp];     // ********** THIS IS MEDIUM

    // DH due to cell adhesion
    // for all neighbours n_i we do Delta( nrg(medium),nrg(n_i)) - Delta(nrg(focal),nrg(n_i) )
    // and sum everything, this amounts to calculate ( Nrg(future) - Nrg(now) )
    for (i = 1; i <= n_nb; i++) {
        int xp2, yp2;
        xp2 = x + nx[i];
        yp2 = y + ny[i];
        if (par.periodic_boundaries) {

            // since we are asynchronic, we cannot just copy the borders once every MCS
            //actually, haven't we dealt with this already in AmoebaeMove? No, this is about all the other neighbours

            if (xp2 <= 0) xp2 = sizex - 2 + xp2;
            if (yp2 <= 0) yp2 = sizey - 2 + yp2;
            if (xp2 >= sizex - 1) xp2 = xp2 - sizex + 2;
            if (yp2 >= sizey - 1) yp2 = yp2 - sizey + 2;

            neighsite = sigma[xp2][yp2];


        } else {
            if (xp2 <= 0 || yp2 <= 0 || xp2 >= sizex - 1 || yp2 >= sizey - 1)
                neighsite = -1;
            else
                neighsite = sigma[xp2][yp2];
        }

        if (neighsite == -1) { // border - for "fixed" boundary conditions
            DH += (sxyp == 0 ? 0 : par.border_energy) - (sxy == 0 ? 0 : par.border_energy);
        } else {
            //DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
            // notice that sxyp is medium, so there is no need of calling the function.
            DH += energyDifference(sxyp, neighsite) - adhesionEnergy(sxy, neighsite);
        }
    }



    // lambda is determined by chemical 0

    //cerr << "[" << lambda << "]";



    //THIS IS THE ONLY CASE HERE -->  if ( sxyp == MEDIUM ) {
    DH += (int) (par.lambda * (1. - 2. * (double) ((*cell)[sxy].Area() - (*cell)[sxy].TargetArea())));
    double ax, ay;
    if ((*cell)[sxy].getMu() > 0.0001 || (*cell)[sxyp].getMu() > 0.0001) {
        double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
        double smeany = (*cell)[sxy].getYpos();

        if (par.periodic_boundaries) {
            if ((x - smeanx) > 0 && (x - smeanx) > (smeanx - (x - (par.sizex - 2)))) {
                smeanx += (par.sizex - 2);
//          cerr<<"dhwm hello"<<endl;
                //cerr<<"passb"<<endl;
            } else if ((smeanx - x) > 0 && (smeanx - x) > (x + (par.sizex - 2) - smeanx)) {
                smeanx -= (par.sizex - 2);
//            cerr<<"dhwm hello"<<endl;
                //cerr<<"passb"<<endl;
            }
            if ((y - smeany) > 0 && (y - smeany) > (smeany - (y - (par.sizey - 2)))) {
                smeany += (par.sizey - 2);
                //cerr<<"passb"<<endl;

            } else if ((smeany - y) > 0 && (smeany - y) > (y + (par.sizey - 2) - smeany)) {
                smeany -= (par.sizey - 2);
                //cerr<<"passb"<<endl;
            }
        }

        ax = x - smeanx;
        ay = y - smeany;
        DH += (*cell)[sxy].getMu() * (ax * (*cell)[sxy].getXvec() + ay * (*cell)[sxy].getYvec()) / hypot(ax, ay);
    }

    //Similarly to Joost's method, a bias due to chemokine gradient
    if ((*cell)[sxy].getChemMu() > 0.0001 || (*cell)[sxyp].getChemMu() > 0.0001) {
        double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
        double smeany = (*cell)[sxy].getYpos();

        if (par.periodic_boundaries) {
            if ((x - smeanx) > 0 && (x - smeanx) > (smeanx - (x - (par.sizex - 2)))) {
                smeanx += (par.sizex - 2);
                //          cerr<<"dhwm hello"<<endl;
                //cerr<<"passb"<<endl;
            } else if ((smeanx - x) > 0 && (smeanx - x) > (x + (par.sizex - 2) - smeanx)) {
                smeanx -= (par.sizex - 2);
                //            cerr<<"dhwm hello"<<endl;
                //cerr<<"passb"<<endl;
            }
            if ((y - smeany) > 0 && (y - smeany) > (smeany - (y - (par.sizey - 2)))) {
                smeany += (par.sizey - 2);
                //cerr<<"passb"<<endl;

            } else if ((smeany - y) > 0 && (smeany - y) > (y + (par.sizey - 2) - smeany)) {
                smeany -= (par.sizey - 2);
                //cerr<<"passb"<<endl;
            }
        }

        ax = x - smeanx;
        ay = y - smeany;
        DH += (*cell)[sxy].getChemMu() * (ax * (*cell)[sxy].getChemXvec() + ay * (*cell)[sxy].getChemYvec()) /
              hypot(ax, ay);


        // cout << "Migrating1!"<<endl;
        //cout<< sxy<<" "<<(*cell)[sxy].getMu()<<" "<<sxyp<<" "<<(*cell)[sxyp].getMu()<<endl;

        // WAS THIS BELOW, BUGGY WITH WRAPPED BOUNDARIES
        //ax=x-(*cell)[sxy].getXpos();
        //ay=y-(*cell)[sxy].getYpos();
        //DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);
    }
    return int(DH);
}


int CellularPotts::DeltaH(int x, int y, int xp, int yp, PDE *PDEfield) {
    double DH = 0.;
    int i, sxy, sxyp;
    int neighsite;

    /* Compute energydifference *IF* the copying were to occur */
    sxy = sigma[x][y];
    sxyp = sigma[xp][yp];

    // cerr<< "x ,y : " << x<<" "<<y<<endl;
    // cerr<< "xp,yp: " << xp<<" "<<yp<<endl;
    //
    /* DH due to cell adhesion */
    for (i = 1; i <= n_nb; i++) {
        int xp2, yp2;
        xp2 = x + nx[i];
        yp2 = y + ny[i];
        if (par.periodic_boundaries) {
            // since we are asynchronic, we cannot just copy the borders once
            // every MCS
            //actually, haven't we dealt with this already in AmoebaeMove? No, this is about all the other neighbours
            if (xp2 <= 0) xp2 = sizex - 2 + xp2;
            if (yp2 <= 0) yp2 = sizey - 2 + yp2;
            if (xp2 >= sizex - 1) xp2 = xp2 - sizex + 2;
            if (yp2 >= sizey - 1) yp2 = yp2 - sizey + 2;
            neighsite = sigma[xp2][yp2];


        } else {
            if (xp2 <= 0 || yp2 <= 0 || xp2 >= sizex - 1 || yp2 >= sizey - 1)
                neighsite = -1;
            else
                neighsite = sigma[xp2][yp2];
        }

        if (neighsite == -1) { // border - for "fixed" boundary conditions
            DH += (sxyp == 0 ? 0 : par.border_energy) - (sxy == 0 ? 0 : par.border_energy);
        } else {
            //cerr<<"Hello DH 0.1"<<endl;

            // MAKE A FUNCTION that takes all values: J's and regulation and sigmas
            //get to the point where you write:
            //DH += adhesionEnergy( s1 , s2, J value)
            //    - adhesionEnergy( same but for the other pair)
            //DH += adhesionEnergy( (*cell)[neighsite].GetExtProtExpress_Fraction() , (*cell)[sxyp].GetExtProtExpress_Fraction() , sxyp , neighsite, (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) )
            //    - adhesionEnergy( (*cell)[neighsite].GetExtProtExpress_Fraction() , (*cell)[sxy].GetExtProtExpress_Fraction() , sxy , neighsite, (*cell)[sxy].EnergyDifference((*cell)[neighsite]) )
            // cerr<< "sigmas sxyp, neigh: "<<sxyp<<" "<<neighsite<<": "<< adhesionEnergy(sxyp , neighsite)<<endl;
            // cerr<< "sigmas sxy, neigh: "<<sxy<<" "<<neighsite<<": "<< adhesionEnergy(sxy , neighsite)<<endl;
            // cerr<< "DH += " << (int) (Adhesion_Energy(sxyp , neighsite) - adhesionEnergy(sxy , neighsite))<<endl ;
            DH += adhesionEnergy(sxyp, neighsite)
                  - adhesionEnergy(sxy, neighsite);
            // NOTICE THAT DH is an integer... dammit! This can create problems when multiplying with a factor
            //DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite])
            //      - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
        }
    }

    // cerr<<"End of Adh, DH = "<< DH<<endl;

    // lambda is determined by chemical 0

    //cerr << "[" << lambda << "]";
    if (sxyp == MEDIUM) {
        DH += (int) (par.lambda * (1. - 2. * (double) ((*cell)[sxy].Area() - (*cell)[sxy].TargetArea())));
    } else if (sxy == MEDIUM) {
        DH += (int) ((par.lambda * (1. + 2. * (double) ((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()))));
    } else
        DH += (int) ((par.lambda * (2. + 2. *
                                         (double) ((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea() -
                                                   (*cell)[sxy].Area() +
                                                   (*cell)[sxy].TargetArea()))));

    //cell migration
    //Joost's method
    double ax, ay;
    if ((*cell)[sxy].getMu() > 0.0001 || (*cell)[sxyp].getMu() > 0.0001) {
        if (sxy != MEDIUM) {
            //cerr<<"tvecx: "<<(*cell)[sxy].getXvec()<<", tvecy: "<< (*cell)[sxy].getYvec() <<endl;
            double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
            double smeany = (*cell)[sxy].getYpos();

            if (par.periodic_boundaries) {
                // if x is on the right and meanx is on the left
                // and if by moving meanx to the right we diminish this distance
                if ((x - smeanx) > 0 && (x - smeanx) > (smeanx + (par.sizex - 2) - x)) {
                    smeanx += (par.sizex - 2);
                    //cerr<<"dh s hello"<<endl;
                    //cerr<<"passb"<<endl;
                } else if ((smeanx - x) > 0 && (smeanx - x) > (x + (par.sizex - 2) - smeanx)) {
                    smeanx -= (par.sizex - 2);
//            cerr<<"dh s hello"<<endl;
                    //cerr<<"passb"<<endl;
                }
                if ((y - smeany) > 0 && (y - smeany) > (smeany - (y - (par.sizey - 2)))) {
                    smeany += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if ((smeany - y) > 0 && (smeany - y) > (y + (par.sizey - 2) - smeany)) {
                    smeany -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }
            }

            ax = x - smeanx;
            ay = y - smeany;
            DH += (*cell)[sxy].getMu() * (ax * (*cell)[sxy].getXvec() + ay * (*cell)[sxy].getYvec()) / hypot(ax, ay);
        }
        if (sxyp != MEDIUM) {
            double spmeanx = (*cell)[sxyp].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
            double spmeany = (*cell)[sxyp].getYpos();

            if (par.periodic_boundaries) {
                if ((x - spmeanx) > 0 && (x - spmeanx) > (spmeanx - (x - (par.sizex - 2)))) {
                    spmeanx += (par.sizex - 2);
//            cerr<<"dh sp hello"<<endl;
                    //cerr<<"passb"<<endl;
                } else if ((spmeanx - x) > 0 && (spmeanx - x) > (x + (par.sizex - 2) - spmeanx)) {
                    spmeanx -= (par.sizex - 2);
//            cerr<<"dh sp hello"<<endl;
                    //cerr<<"passb"<<endl;
                }
                if ((y - spmeany) > 0 && (y - spmeany) > (spmeany - (y - (par.sizey - 2)))) {
                    spmeany += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if ((spmeany - y) > 0 && (spmeany - y) > (y + (par.sizey - 2) - spmeany)) {
                    spmeany -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }
            }
            //ax=x-(*cell)[sxyp].getXpos(); //returns meanx
            //ay=y-(*cell)[sxyp].getYpos(); //returns meany
            ax = x - spmeanx;
            ay = y - spmeany;
            DH -= (*cell)[sxyp].getMu() * (ax * (*cell)[sxyp].getXvec() + ay * (*cell)[sxyp].getYvec()) / hypot(ax, ay);
        }
    }


    //Similar to Joost's method, but for chemotaxis (no persistence!)
    if ((*cell)[sxy].getChemMu() > 0.0001 || (*cell)[sxyp].getChemMu() > 0.0001) {
        if (sxy != MEDIUM) {
            //cerr<<"tvecx: "<<(*cell)[sxy].getXvec()<<", tvecy: "<< (*cell)[sxy].getYvec() <<endl;
            double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
            double smeany = (*cell)[sxy].getYpos();

            if (par.periodic_boundaries) {
                // if x is on the right and meanx is on the left
                // and if by moving meanx to the right we diminish this distance
                if ((x - smeanx) > 0 && (x - smeanx) > (smeanx + (par.sizex - 2) - x)) {
                    smeanx += (par.sizex - 2);
                    //cerr<<"dh s hello"<<endl;
                    //cerr<<"passb"<<endl;
                } else if ((smeanx - x) > 0 && (smeanx - x) > (x + (par.sizex - 2) - smeanx)) {
                    smeanx -= (par.sizex - 2);
                    //            cerr<<"dh s hello"<<endl;
                    //cerr<<"passb"<<endl;
                }
                if ((y - smeany) > 0 && (y - smeany) > (smeany - (y - (par.sizey - 2)))) {
                    smeany += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if ((smeany - y) > 0 && (smeany - y) > (y + (par.sizey - 2) - smeany)) {
                    smeany -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }
            }

            ax = x - smeanx;
            ay = y - smeany;
            DH += (*cell)[sxy].getChemMu() * (ax * (*cell)[sxy].getChemXvec() + ay * (*cell)[sxy].getChemYvec()) /
                  hypot(ax, ay);
        }
        if (sxyp != MEDIUM) {
            double spmeanx = (*cell)[sxyp].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other length
            double spmeany = (*cell)[sxyp].getYpos();

            if (par.periodic_boundaries) {
                if ((x - spmeanx) > 0 && (x - spmeanx) > (spmeanx - (x - (par.sizex - 2)))) {
                    spmeanx += (par.sizex - 2);
                    //            cerr<<"dh sp hello"<<endl;
                    //cerr<<"passb"<<endl;
                } else if ((spmeanx - x) > 0 && (spmeanx - x) > (x + (par.sizex - 2) - spmeanx)) {
                    spmeanx -= (par.sizex - 2);
                    //            cerr<<"dh sp hello"<<endl;
                    //cerr<<"passb"<<endl;
                }
                if ((y - spmeany) > 0 && (y - spmeany) > (spmeany - (y - (par.sizey - 2)))) {
                    spmeany += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if ((spmeany - y) > 0 && (spmeany - y) > (y + (par.sizey - 2) - spmeany)) {
                    spmeany -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }
            }
            ax = x - spmeanx;
            ay = y - spmeany;
            DH -= (*cell)[sxyp].getChemMu() * (ax * (*cell)[sxyp].getChemXvec() + ay * (*cell)[sxyp].getChemYvec()) /
                  hypot(ax, ay);
        }
    }
    // This could be a double but some optimizations in CopyProb require it to be an int
    return int(round(DH));
}


double CellularPotts::energyDifference(int sigma1, int sigma2) const {
    if (not (sigma1 and sigma2)) {
        return par.Jmed;
    }
    const Cell &c1 = (*cell)[sigma1];
    const Cell &c2 = (*cell)[sigma2];
    if (c1.getGroup() != c2.getGroup())
        return par.Jmed * 2;
    return par.Jalpha
           + getKLStrength(c1.jkey_dec, c2.jlock_dec)
           + getKLStrength(c2.jkey_dec, c1.jlock_dec);
}


double CellularPotts::adhesionEnergy(int sigma1, int sigma2) const {
    if (sigma1 == sigma2) return 0;
    return energyDifference(sigma1, sigma2);
}


double CellularPotts::calculateGamma(int sigma1, int sigma2) const {
    // energyDifference is used instead of AdhesionStrength because we might want to calculate c1 x c1
    return par.Jmed - energyDifference(sigma1, sigma2) / 2;
}


//this the function that changes the getSigma
// it also does book-keeping of everything
void CellularPotts::ConvertSpinToMedium(int x, int y) {
    int tmpcell;
    if ((tmpcell = sigma[x][y])) { // if tmpcell is not MEDIUM - this should be excluded outside, so it is good check
        (*cell)[tmpcell].DecrementArea(); // returns --area
        (*cell)[tmpcell].RemoveSiteFromMoments(x, y);

        //if cell's area is zero -> kill it
        if (!(*cell)[tmpcell].Area()) {
            (*cell)[tmpcell].Apoptose();
            //cerr << "Cell " << tmpcell << " apoptosed\n";
        }
    } else {
        cerr << "ConvertSpinToMedium(): Error. Trying to convert spin of medium with medium" << endl;
    }

//   if ( (tmpcell=getSigma[xp][yp]) ) {// if tmpcell is not MEDIUM, it gains one pixel
//     (*cell)[tmpcell].IncrementArea();
//     (*cell)[tmpcell].AddSiteToMoments(x,y);
//
//   }

    // Book-keeping of contacts
    int check = 0;
    int point;
    int idSelf = sigma[x][y];     // Surely not medium <- we check this in this special case
    int idNeigh = MEDIUM;  // MEDIUM IN MY CASE

    ///use this if you need info about contacts with neighbours
    // go through neighbourhood
    for (int k = 1; k <= n_nb; k++) {
        // Do you not update contacts at boundaries???
        int neix = x + nx[k];
        int neiy = y + ny[k];
        if (neix <= 0 || neix >= sizex - 1 || neiy <= 0 || neiy >= sizey - 1) {
            if (par.periodic_boundaries) {
                if (neix <= 0) neix = sizex - 2 + neix;
                if (neix >= sizex - 1) neix = neix - sizex + 2;
                if (neiy <= 0) neiy = sizey - 2 + neiy;
                if (neiy >= sizey - 1) neiy = neiy - sizey + 2;
            } else {
                continue;
            }
        }


        //sigmaneigh = getSigma[ neix ][ neiy ];



        //if(x+nx[k]>0 && x+nx[k]<sizex-1 && y+ny[k]>0 && y+ny[k]<sizey-1 ){
        //  point=sigma[x+nx[k]][y+ny[k]]; //take getSigma of neighbour, this can be medium
        point = sigma[neix][neiy]; //take getSigma of neighbour, this can be medium
        // if neigh is not the same as self (you don't have boundaries with yourself, are you thereby really free alone? :P )
        // if neigh is not same as self, there are a few cases:
        // self is not medium <- we update
        if (point != idSelf) {
            //if idSelf not medium, we will remove this point from its contact (because we are going to copy neigh in it)
            if (idSelf)
                check += (*cell)[idSelf].updateNeighbourBoundary(point, -1);
            if (point)
                check += (*cell)[point].updateNeighbourBoundary(idSelf,
                                                                -1); //and if point is not medium, we are removing self from contacts of neigh
        }
        // if neigh is not the same as
        if (point != idNeigh) {
            if (idNeigh)
                check += (*cell)[idNeigh].updateNeighbourBoundary(point, 1);
            if (point)
                check += (*cell)[point].updateNeighbourBoundary(idNeigh, 1);
        }
        if (check) {
            printf(
                "error in ConvertSpinToMedium(): wrongly updating neighbours of copy event idSelf %d and idNeigh %d\n",
                idSelf, idNeigh);
            //printf("agent nr %d\n", agentid);
            exit(1);
        }
        //}
    }

    //cout<<"in ConvertSpin, after: "<<(*cell)[1].neighbours[0].first<<endl;


    sigma[x][y] = MEDIUM;    // THE SPIN COPYING
    //cerr<<"FYI: ConvertSpinToMedium() happened" << endl;
    //exit(1);

}


//this the function that changes the getSigma
// it also does book-keeping of everything
void CellularPotts::ConvertSpin(int x, int y, int xp, int yp) {
    int tmpcell;
    if ((tmpcell = sigma[x][y])) { // if tmpcell is not MEDIUM
        (*cell)[tmpcell].DecrementArea(); // returns --area
        (*cell)[tmpcell].RemoveSiteFromMoments(x, y);

        //if cell's area is zero -> kill it
        if (!(*cell)[tmpcell].Area()) {
            (*cell)[tmpcell].Apoptose();
            //cerr << "Cell " << tmpcell << " apoptosed\n";
        }
    }

    if ((tmpcell = sigma[xp][yp])) {// if tmpcell is not MEDIUM, it gains one pixel
        (*cell)[tmpcell].IncrementArea();
        (*cell)[tmpcell].AddSiteToMoments(x, y);

    }

    // Book-keeping of contacts
    int check = 0;
    int point;
    int idSelf = sigma[x][y];
    int idNeigh = sigma[xp][yp];
    ///use this if you need info about contacts with neighbours
    for (int k = 1; k <= n_nb; k++) {
        // Do you not update contacts at boundaries???
        int neix = x + nx[k];
        int neiy = y + ny[k];
        if (neix <= 0 || neix >= sizex - 1 || neiy <= 0 || neiy >= sizey - 1) {
            if (par.periodic_boundaries) {
                if (neix <= 0) neix = sizex - 2 + neix;
                if (neix >= sizex - 1) neix = neix - sizex + 2;
                if (neiy <= 0) neiy = sizey - 2 + neiy;
                if (neiy >= sizey - 1) neiy = neiy - sizey + 2;
            } else {
                continue;
            }
        }


        //sigmaneigh = getSigma[ neix ][ neiy ];



        //if(x+nx[k]>0 && x+nx[k]<sizex-1 && y+ny[k]>0 && y+ny[k]<sizey-1 ){
        //  point=sigma[x+nx[k]][y+ny[k]]; //take getSigma of neighbour, this can be medium
        point = sigma[neix][neiy]; //take getSigma of neighbour, this can be medium

        //Error goes that point = 658 = idNeigh, idSelf=MEDIUM
        if (point != idSelf) {
            if (idSelf) {
                check += (*cell)[idSelf].updateNeighbourBoundary(point, -1);
                if (check) cerr << "Here1" << endl;
            }
            if (point) {
                check += (*cell)[point].updateNeighbourBoundary(idSelf, -1);
                if (check) cerr << "Here2" << endl;
            }
        }
        if (point != idNeigh) {
            if (idNeigh)
                check += (*cell)[idNeigh].updateNeighbourBoundary(point, 1);
            if (point)
                check += (*cell)[point].updateNeighbourBoundary(idNeigh, 1);
        }
        if (check) {
            printf("error in ConvertSpin(): wrongly updating neighbours of copy event idSelf %d and idNeigh %d\n",
                   idSelf,
                   idNeigh);
            cout << "Extra info: idSelf=sigma[x][y], where x and y are " << x << "," << y << endl;
            cout << "Extra info: idNeigh=sigma[xp][yp], where xp and yp are " << xp << "," << yp << endl;
            cout << "Extra info: neix,neiy " << neix << "," << neiy << endl;
            //printf("agent nr %d\n", agentid);

            exit(1);
        }
        //}
    }

    //cout<<"in ConvertSpin, after: "<<(*cell)[1].neighbours[0].first<<endl;

    sigma[x][y] = sigma[xp][yp]; // = idNeigh - THE SPIN COPYING


}


/** PUBLIC **/
// stiff = conn_diss i.e. the Delta H coming from break connectivity of cell, it is zero in debugging
int CellularPotts::CopyvProb(int DH, double stiff) {

    double dd;
    int s;
    s = (int) stiff;
    if (DH <= -s) return 2; //accept copy if DH+s <= 0

    // if DH becomes extremely large, calculate probability on-the-fly
    // common values - from 0 to 1023 - are stored in copyprob - checked
    // BOLTZMANN = 1024 and is defined in sticky.h
    if (DH + s > BOLTZMANN - 1)
        dd = exp(-((double) (DH + s) / par.T));
    else
        dd = copyprob[DH + s];

    //for(int i=0;i<1024;i++) cerr<<i<<" "<<copyprob[i]<<endl;
    //exit(1);

    if (RANDOM() < dd) return 1; else return 0;
}

void CellularPotts::CopyProb(double T) {
    int i;
    for (i = 0; i < BOLTZMANN; i++)
        copyprob[i] = exp(-((double) (i) / T));
}

#include <fstream>

//! Monte Carlo Step. Returns summed energy change
// In this version there is a chance that you copy medium (from outer space)
int CellularPotts::AmoebaeMove2(PDE *PDEfield) {
    double chanceofmedium = par.chancemediumcopied;  //this is the chance that we test for medium instead of actual pixel
    // (i.e. that we raise mediumflag for that)
    //this is NOT probability that medium is inserted,
    // that follows the normal process and depends on nrg.
    //BTW, we raise this after we are sure that k != kp: we don't want medium inside the cell :P

    bool mediumflag, converted, inout1, inout2;

    double loop, p;
    //int updated=0;
    thetime++;
    int SumDH = 0;

    if (frozen)
        return 0;

    //variables used in loop
    int x, y, xp, yp, xn, yn;
    int k, kp;
    int edgesize, rindex;

    loop = (double) edgeSetVector.size_map() / n_nb;
    //cout <<" loop size "<<loop<<endl;
    for (int i = 0; i < loop; i++) {
        mediumflag = false;   // what happens if this is not re-set to false all the times?
        // that the first time it goes true it will be true for the rest of the loop!
        converted = false;
        edgesize = edgeSetVector.size_map(); //are there still edges to update?
        if (edgesize == 0) { break; }
        else {
            rindex = RandomNumber(edgesize) - 1; //pick random edge (Randomnumber goes from 1-N)
            auto it = edgeSetVector[rindex]; //this returns an iterator to an unordered map element out of the vector edgevector in
            x = (*it).first[0];
            y = (*it).first[1];
            xp = (*it).first[2];
            yp = (*it).first[3];
            k = sigma[x][y];

            if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
                kp = -1;
            else
                kp = sigma[xp][yp];
            // test for border state (relevant only if we do not use periodic boundaries)
            // test always passed with periodic boundaries
            if (kp != -1) {
                // Don't even think of copying the special border state into you!
                if (k != kp) {
                    // connectivity dissipation:
                    int H_diss = 0;

                    if (!ConnectivityPreservedP(x, y))
                        H_diss = par.conn_diss;

                    //if k is different from medium and with small chance = chancemedium,
                    // we propose to copy medium into k
                    if (k != MEDIUM && RANDOM() < chanceofmedium) {
                        kp = MEDIUM; //and so kp = MEDIUM
                        mediumflag = true;
                    }
                    //IF we are not copying medium from outer space
                    if (!mediumflag) {
                        int D_H = DeltaH(x, y, xp, yp, PDEfield);

                        if ((p = CopyvProb(D_H, H_diss)) > 0) {
                            ConvertSpin(x, y, xp, yp);
                            converted = true;
                            SumDH += D_H;
                        }
                    } else {
                        int D_H = DeltaHWithMedium(x, y, PDEfield);

                        if ((p = CopyvProb(D_H, H_diss)) > 0) {
                            ConvertSpinToMedium(x, y);
                            converted = true;
                            SumDH += D_H;
                        }

                    }
                    if (converted) {
                        //go through neighbourhood
                        for (int j = 1; j <= n_nb; j++) {
                            xn = nx[j] + x;
                            yn = ny[j] + y;

                            if (xn > 0 && yn > 0 && xn < sizex - 1 &&
                                yn < sizey - 1) {//if the neighbour site is within the lattice
                                auto kn = sigma[xn][yn];
                                // if we should add the edge to the edgelist, add it
                                if (kn == kp) { //if we should add the edge to the edgelist, add it
                                    edgeSetVector.erase({x, y, xn, yn});
                                    edgeSetVector.erase({xn, yn, x, y});
                                    loop -= 2. / n_nb;
                                } else if (kn >= 0) {
                                    inout1 = edgeSetVector.insert({x, y, xn, yn});
                                    inout2 = edgeSetVector.insert({xn, yn, x, y});
                                    loop += (inout1 + inout2) / (double) n_nb;
                                }
                            } //end if neighbour in lattice
                        }//end neighbourhood loop
                    }//end if converted

                } //end if neighbour is different pixel type
            }//end if neighbour is in the field
        } //if edgelist still full
    }//end for loop

    return SumDH;

}


// void CellularPotts::ReadZygotePicture(void) {
//
//
//
//   int pix,cells,i,j,c,p,checkx,checky;
//   char **pixelmap;
//   char pixel[3];
//
//   sscanf(ZYGXPM(ZYGOTE)[0],"%d %d %d %d",&checkx,&checky,&cells,&pix);
//
//   if ((checkx>sizex)||(checky>sizey)) {
//     std::cerr <<  "ReadZygote: The included xpm picture is smaller than the grid!\n";
//     std::cerr << "\n Please adjust either the grid size or the picture size.\n";
//     std::cerr << sizex << "," << sizey << "," << checkx << "," << checky << "\n";
//     exit(1);
//   }
//
//   pixelmap=(char **)malloc(cells*sizeof(char *));
//   if (pixelmap==NULL) MemoryWarning();
//
//   pixelmap[0]=(char *)malloc(cells*3*sizeof(char));
//   if (pixelmap[0]==NULL) MemoryWarning();
//
//   for(i=1;i<cells;i++)
//     pixelmap[i]=pixelmap[i-1]+3;
//
//   for (i=0;i<cells;i++) {
//     for (j=0;j<pix;j++)
//       pixelmap[i][j]=ZYGXPM(ZYGOTE)[i+1][j];
//     pixelmap[i][pix]='\0';
//   }
//
//   for (i=0;i<sizex*sizey;i++) getSigma[0][i]=0;
//   fprintf(stderr,"[%d %d]\n",checkx,checky);
//
//   int offs_x, offs_y;
//   offs_x=(sizex-checkx)/2;
//   offs_y=(sizey-checky)/2;
//
//   for (i=0;i<checkx;i++)
//     for (j=0;j<checky;j++) {
//       for (p=0;p<pix;p++)
//         pixel[p]=ZYGXPM(ZYGOTE)[cells+1+j][i*pix+p];
//
//       pixel[pix]='\0';
//
//       for (c=0;c<cells;c++) {
// 	if (!(strcmp(pixelmap[c],pixel))) {
// 	  if ( (getSigma[offs_x+i][offs_y+j]=c) ) {
//
// 	    // if c is _NOT_ medium (then c=0)
// 	    // assign pixel values from "sigmamax"
// 	    getSigma[offs_x+i][offs_y+j]+=(Cell::MaxSigma()-1);
// 	  }
// 	}
//
//       }
//     }
//
//   free(pixelmap[0]);
//   free(pixelmap);
// }


void CellularPotts::ConstructInitCells(Dish &beast) {

    // Get the maximum cell ID (mostly equal to the cell number)
    //int loop=sizex*sizey;
    int cells = 0;
    for (int i = 1; i < sizex - 1; i++)
        for (int j = 1; j < sizey - 1; j++) {
            if (cells < sigma[i][j]) cells = sigma[i][j];
        }
    cout << "nr cells placed " << cells << endl;

    cout << "[ cells = " << cells << "]\n";

    // construct enough cells for the zygote.  "cells", contains the
    // number of colours (excluding background).
    {
        for (int i = 0; i < cells; i++) {
            cell->push_back(Cell(beast));
        }
    }

    //initialises mean x and y to some values, does not matter what,
    // as long as meanx and y get a defined number BEFORE we calculate MeasureCellSize()
    // which, in turn, if periodic_boundaries = true, depends on meanx and y
    for (auto &c: *cell) {
        c.InitMeanX(sizex / 2.);
        c.InitMeanY(sizey / 2.);
    }


    // Set the area and target area of the cell
    // makes use of the pointer to the Cell pointer of Dish
    // which is a member of CellularPotts
    MeasureCellSizes();

    //for (vector<Cell>::iterator c=cell->begin(); c!=cell->end();c++) {
    //  cerr<<"getSigma: "<<c->Sigma()<<". New meanx and y calculated as: "<<c->getXpos()<<" "<<c->getYpos()<<endl;
    //}

    // set zygote_area to mean cell area.
    int mean_area = 0;
    for (auto &c: *cell) {
        mean_area += c.Area();
    }
    if (cells != 0)
        mean_area /= cells;

    zygote_area = mean_area;

    cout << "mean_area = " << mean_area << "\n";
    // set all cell areas to the mean area
    {
        for (auto &c: *cell) {
            c.setMu(0.0);

            if (par.target_area) {
                c.SetTargetArea(par.target_area);
            } else {
                c.SetTargetArea(mean_area);

            }
        }
    }
    cout << "ConstructInitCells is done" << endl;
}

void CellularPotts::MeasureCellSizes() {

    // Clean areas of all cells, including medium
    for (auto &c: *cell) {
        c.SetTargetArea(0);
        c.area = 0;
        c.CleanMoments();
    }

    // calculate the area of the cells
    for (int x = 1; x < sizex - 1; x++) {
        for (int y = 1; y < sizey - 1; y++) {
            if (sigma[x][y]) {
                (*cell)[sigma[x][y]].IncrementTargetArea();
                (*cell)[sigma[x][y]].IncrementArea();
                (*cell)[sigma[x][y]].AddSiteToMoments(x, y);

            }
        }
    }

    // set the actual area to the target area
    {
        for (auto &c: *cell) {
            c.SetAreaToTarget();
        }
    }
}

void CellularPotts::MeasureCellSize(Cell &c) {

    c.CleanMoments();

    // calculate the area of the cell
    for (int x = 1; x < sizex - 1; x++) {
        for (int y = 1; y < sizey - 1; y++) {
            if (sigma[x][y] == c.sigma) {
                (*cell)[sigma[x][y]].IncrementTargetArea();
                (*cell)[sigma[x][y]].IncrementArea();
                (*cell)[sigma[x][y]].AddSiteToMoments(x, y);

            }
        }
    }

//   // set the actual area to the target area
//   {
//   for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
//     c->SetAreaToTarget();

//   }

}

// Rewritten version of the funciton below
// it uses the cell meanx and meany for the center of mass of the cell,
// IN THE FUTURE -> only calculates PCA for cells that need to be split...
Dir *CellularPotts::FindCellDirections3() const {
    // double array - allocated (why c style?) the size population
    double *sumx = nullptr, *sumy = nullptr;
    double *sumxx = nullptr, *sumxy = nullptr, *sumyy = nullptr;
    double *n = nullptr;

    double xmean = 0, ymean = 0, sxx = 0, sxy = 0, syy = 0;
    double D, lb1 = 0, lb2 = 0;

    Dir *celldir;


    /* Allocation of sufficient memory space */
    if ((sumx = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();
    if ((sumy = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();
    if ((sumxx = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();
    if ((sumxy = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();
    if ((sumyy = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();

    if ((n = (double *) malloc(cell->size() * sizeof(double))) == nullptr)
        MemoryWarning();

    if (!(celldir = new Dir[cell->size()]))
        MemoryWarning();


    /* Initialization of the variables */
    for (int i = 0; i < (int) cell->size(); i++) {
        sumx[i] = 0.;
        sumy[i] = 0.;
        sumxx[i] = 0.;
        sumxy[i] = 0.;
        sumyy[i] = 0.;
        n[i] = 0L;
    }

    /* Find sumx, sumy, sumxx and sumxy for all cells */
    for (int x = 1; x < sizex - 1; x++) {
        for (int y = 1; y < sizey - 1; y++) {
            //for (int x=0;x<sizex;x++)
            //  for (int y=0;y<sizey;y++)
            if (sigma[x][y] > 0) {
                //cerr<<getSigma[x][y]<<endl;

                sumx[0] += (double) x;  //pos 0 contains global statistics
                sumy[0] += (double) y;
                sumxx[0] += (double) x * x;
                sumxy[0] += (double) x * y;
                sumyy[0] += (double) y * y;
                n[0]++;

                int isigma = sigma[x][y];

                double tmpx = x; //because we may change them if wrapped
                double tmpy = y;

                double meanx = (*cell)[isigma].meanx;
                double meany = (*cell)[isigma].meany;

                if (par.periodic_boundaries) {
                    //Now the wrapping - in this function we wrap around meanx:
                    //if x is closer to running average when wrapped, we wrap it
                    if (x - meanx > 0 && x - meanx > (meanx + (sizex - 2) - x)) {
                        tmpx -= (sizex - 2);
                        //cerr<<"directions passb1"<<endl;
                    } else if (meanx - x > 0 && meanx - x > (x + (sizex - 2) - meanx)) {
                        tmpx += (sizex - 2);
                        //cerr<<"directions passb1"<<endl;
                    }
                    //same for y
                    if (y - meany > 0 && y - meany > (meany + (sizey - 2) - y)) {
                        tmpy -= (sizey - 2);
                        //cerr<<"directions passb1"<<endl;
                    } else if (meany - y > 0 && meany - y > (y + (sizey - 2) - meany)) {
                        tmpy += (sizey - 2);
                        //cerr<<"directions passb1"<<endl;
                    }
                }

                sumx[isigma] += tmpx; //pos indicised by getSigma - local stats about x direction
                sumy[isigma] += tmpy;

                //Different model here, we just do <(x-meanx)^2>
                sumxx[isigma] += (pow(tmpx - meanx, 2.));
                sumyy[isigma] += (pow(tmpy - meany, 2.));
                sumxy[isigma] += ((tmpx - meanx) * (tmpy - meany));

                // This is how it was
                //sumxx[isigma]+=(double)tmpx*tmpx;
                //sumxy[isigma]+=(double)tmpx*tmpy; // mixed coordinate
                //sumyy[isigma]+=(double)tmpy*tmpy;

                n[isigma]++;

            }
        }
    }
    //exit(1);
    // Compute the principal axes for all cells from eigenvalues of correlation matrix C
    //     ( sxx sxy )
    // C = <          >
    //     ( sxy syy )
    // We diagonalise this and find eigenvalues lb1 and lb2
    //recalculate the means while we're at it


    double small_number = 0.0000001;
    double large_enough_number = 1. / small_number;

    for (int i = 0; i < (int) cell->size(); i++) {
        celldir[i].meanx = (*cell)[i].meanx;
        celldir[i].meany = (*cell)[i].meany;
        //cerr<<"celldir[i].bb1 = "<<celldir[i].bb1<<endl;  // there are some problems here and at initilisation (see video)
        double xmean = celldir[i].meanx;
        double ymean = celldir[i].meany;
        sxx = 0.;
        syy = 0.;
        sxy = 0.;
        // why this? of course moments are ill defined for small cells
        // maybe it's no problem because they don't divide

        if (n[i] > 10) {
            sxx = sumxx[i] / ((double) n[i]); //or maybe n[i]-1 ? Would be strange because this is all the data
            syy = sumyy[i] / ((double) n[i]);
            sxy = sumxy[i] / ((double) n[i]);

            // Now we have the elements of the matrix, we diagonalise by solving eigenvalue problem (call x the eigenvalues)
            // x^2 - (sxx+syy)x + (sxx*syy-sxy*sxy) = 0
            D = sqrt((sxx + syy) * (sxx + syy) - 4. * (sxx * syy - sxy * sxy)); // this is the discriminant
            lb1 = (sxx + syy + D) / 2.;    // these are the two solutions, i.e. the two eigenvalues
            lb2 = (sxx + syy - D) / 2.;
            // Now lb1 > lb2 because ... well, look at it, then lb1 is largest eigenvalue,
            // so its eigenvector is the principal axis of the cell
            celldir[i].lb1 = lb1;
            celldir[i].lb2 = lb2;
        }
        //Now we get the eigenvectors:
        // first eigenvector v1 is the solution of C - lb1*v1 =0
        // C is covar. matrix, lb1 is eigenv 1 just calculated.
        // v1 has an x and a y component, here ve express v1 as x/y = sxy/(lb1-sxx)

        //this case is when there is no covariance, so the cartesian axis are already the best basis
        if (sxy != 0.0) {
//       cerr<<"sxy!=0, lb1: "<<lb1<<", syy: "<<syy<<", sxx: "<<sxx<<", sxy: "<<sxy<<endl;
            celldir[i].bb1 =
                sxy / (lb1 - syy); // APPARENTLY THIS IS FINE- This is y/x - shoudn't it be = syy/(sxy-lb1)?
            // I think this is correct : celldir[i].bb1=syy/(lb1-sxy);
            // and NOT WHAT IS THERE.. UNLESS I am wrong :P
            if (fabs(celldir[i].bb1) < small_number) {
                if (celldir[i].bb1 > 0.)
                    celldir[i].bb1 = small_number;
                else
                    celldir[i].bb1 = -small_number;
            }

            celldir[i].aa1 =
                ymean -
                xmean * celldir[i].bb1; //this is the intercept to the y axis of the line with slope first eigenv.
            // which passes through xmean and ymean
            celldir[i].bb2 = (-1.) / celldir[i].bb1; // bb2 is the direction perpendicular to the first eigenvector
            // (because the perpend. to a line y=mx+q has slope -1/m)

            celldir[i].aa2 =
                ymean -
                celldir[i].bb2 * xmean; // this is the intercept to y axis of the line of the second eigenvector
            // along this line we cut the cell !!!
        } else {
//       cerr<<"sxy=0, ";
            // USED TO BE
            //celldir[i].bb1=1.; WHICH IS DEFINITELY WRONG

            //Because later we are doing operations on the slope, we should choose a large enough number that does not overflow
            // a good idea could be to choose a slope so that the almost vertical line makes less than epsilon=0.1 error across the whole field
            // so that division is effectively vertical (or horizontal)
            // this is a line that has slope par.sizex/epsilon, let's add a little bit to be extra safe (times 2)
            // with a field size of 1000 and epsilon = 0.1 -> m= 10000 that's ok small for more calculations

            //double large_enough_number = (2.*(double)par.sizex)/0.1;
            //double large_enough_number = 1./0.0000001;
            double random_plus_or_minus_1 = -1 + 2 * (int) (2. * RANDOM());
            if (sxx > syy) {
                celldir[i].bb1 = 0.;
                celldir[i].aa1 = ymean;
                celldir[i].bb2 = random_plus_or_minus_1 * large_enough_number;
                celldir[i].aa2 = ymean - celldir[i].bb2 * xmean;
//         cerr<<"sxx>syy"<<endl;
            } else if (syy > sxx) {
                celldir[i].bb1 = random_plus_or_minus_1 * large_enough_number;
                celldir[i].aa1 = ymean - xmean * celldir[i].bb1;
                //celldir[i].bb2 = 0.;
                celldir[i].bb2 = -1. * random_plus_or_minus_1 * small_number;
                celldir[i].aa2 = ymean;
//         cerr<<"syy>sxx"<<endl;
            } else {
//         cerr<<"syy=sxx"<<endl;
                celldir[i].bb1 = (RANDOM() < 0.5) ? 0. : (random_plus_or_minus_1 *
                                                          large_enough_number); //if sxx==syy we randomise vertical or horizontal
                if (celldir[i].bb1 > 1. || celldir[i].bb1 < 1.) {
                    celldir[i].aa1 = ymean - xmean * celldir[i].bb1;
                    //celldir[i].bb2 = 0.;
                    celldir[i].bb2 = -1. * random_plus_or_minus_1 * small_number;
                    celldir[i].aa2 = ymean;
                } else {
                    celldir[i].aa1 = ymean;
                    celldir[i].bb2 = random_plus_or_minus_1 * large_enough_number;
                    celldir[i].aa2 = ymean - celldir[i].bb2 * xmean;
                }
            }
        }
//     cerr<<"Sigma: "<<i<<", bb2: "<<celldir[i].bb2<<", aa2: "<<celldir[i].aa2<<endl;
    }

    //}

    /* bevrijd gealloceerd geheugen */
    free(sumx);
    free(sumy);
    free(sumxx);
    free(sumxy);
    free(sumyy);
    free(n);

    return celldir;

}


// TODO: This should probably be in cell.h
int CellularPotts::DivideCell(int cell_sigma, BoundingBox box) {

    int sigmaneigh;

    // for the cell directions
    Dir *celldir = nullptr;
    vector<int> toprint;

    Cell *daughterp;
    Cell *motherp = &((*cell)[cell_sigma]);
    motherp->food /= 2;
    /* division */

    //we first check if we can recycle some position already exisiting in the vector
    //such position would come from a cell that has previously apoptosed
    vector<Cell>::iterator c;
    bool replaced = false;
    for (c = cell->begin(), c++; c != cell->end(); c++) {
        if (!c->AliveP() && c->TargetArea() <= 0 && c->Area() == 0) {
            //we recycle this getSigma for the new cell
            //set recycled getSigma
            daughterp = new Cell(*(motherp->owner), motherp->getTau(), c->Sigma());
            daughterp->CellBirth(*motherp);
            *c = *daughterp;    // notice that the operator = (equal) is overloaded, see cell.h
            replaced = true;
            break;
        }
    }
    if (!replaced) {
        //cout << "We do not recycle, calling function new"<< endl;
        // THIS USED TO BE ABOVE, WHERE THE SIGN *** IS !!!
        //MAKES NEW CELL AT THE END OF ARRAY
        daughterp = new Cell(*(motherp->owner));  //this calls  Cell(...){ plane=&who; ConstructorBody()}
        int momcol = motherp->group;
        daughterp->CellBirth(*motherp);
        cell->push_back(*daughterp);  //this calls default copy constructor Cell(const Cell &src)
        // prints "Tomato"
        //this puts new cells at the end of array if there was no space to recycle
        // renew pointer to mother (because after push_back memory might be relocated)
        motherp = &((*cell)[cell_sigma]);
        if (motherp->group != momcol) cerr << "this is the problem" << endl;
    }
    // renew pointers
    delete daughterp;
    if (replaced) {
        daughterp = &(*c);
        //cerr<<"mother getSigma: "<<motherp->Sigma()<<", daughter sigma"<<daughterp->Sigma()<<endl;
    } else {
        //cerr<<"mother sigma: "<<motherp->Sigma()<<", daughter getSigma"<<daughterp->Sigma()<<endl;
        daughterp = &(cell->back());
    }

    int pixel_count = 0;
    for (int i = box.getMinX(); i <= box.getMaxX(); i++)
        for (int j = box.getMinY(); j <= box.getMaxY(); j++) {
            if (sigma[i][j] == motherp->sigma) {
                ++pixel_count;
                /* Now the actual division takes place */

                /* If celldirections where not yet computed: do it now */
                if (!celldir)
                    celldir = FindCellDirections3();
                // if site is below the minor axis of the cell: getSigma of new cell
                // to properly choose this we have to check where this pixel is
                int checki = i;
                int checkj = j;

                if (checkj > ((int) (celldir[motherp->sigma].aa2 + celldir[motherp->sigma].bb2 * (double) checki))) {
                    motherp->DecrementArea();
                    motherp->RemoveSiteFromMoments(i, j);
                    sigma[i][j] = daughterp->Sigma();
                    daughterp->IncrementArea();
                    daughterp->AddSiteToMoments(i, j);
                    //go through neighbourhood to update contacts
                    // to new daughter contacts we now pass duration from mother
                    // getSigma[i][j] is daughter, sigmaneigh can be daughter, mother, medium, someone else
                    for (int k = 1; k <= n_nb; k++) {
                        //if wrapped boundaries we wrap i+nx[k] and j+ny[k] around (if needed)
                        //if fiexed boundaries, we exclude them from the neigh counting
                        int neix = i + nx[k];
                        int neiy = j + ny[k];
                        sigmaneigh = sigma[neix][neiy];
                        if (sigmaneigh < 0 || !BoundingBox(0, 0, sizex - 1, sizey - 1).inside(neix, neiy))
                            continue;
                        //if sigmaneigh is not getSigma, we update the contact of daughter cell with it,
                        //and the contact of that cell with daughter (provided it is not medium)
                        if (sigmaneigh != sigma[i][j]) {

                            //update the edgeSetVector
                            edgeSetVector.insert({i, j, neix, neiy});
                            edgeSetVector.insert({neix, neiy, i, j});

                            (*cell)[sigma[i][j]].updateNeighbourBoundary(sigmaneigh, 1);
                            //take duration from mother iff sigmaneigh is not mother
                            if (sigmaneigh != motherp->Sigma() && sigmaneigh != MEDIUM)
                                (*cell)[sigma[i][j]].SetNeighbourDurationFromMother(sigmaneigh, motherp->returnDuration(
                                    sigmaneigh));

                            //also cell to which sigmaneigh belongs must be updated, if it is not medium
                            if (sigmaneigh) {
                                (*cell)[sigmaneigh].updateNeighbourBoundary(sigma[i][j], 1);
                                if (sigmaneigh != motherp->Sigma()) {
                                    (*cell)[sigmaneigh].SetNeighbourDurationFromMother(sigma[i][j],
                                                                                       motherp->returnDuration(
                                                                                           sigmaneigh));
                                }
                            }

                            if (sigmaneigh != motherp->Sigma()) {
                                motherp->updateNeighbourBoundary(sigmaneigh, -1);

                                if (sigmaneigh)
                                    (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(), -1);
                            }
                        } else//sigmaneigh==getSigma[i][j] This pixel has already become a daughter pixel,
                            //remove from contacts between mother and daughter
                        {
                            //update the edgeSetVector
                            edgeSetVector.erase({i, j, neix, neiy});
                            edgeSetVector.erase({neix, neiy, i, j});

                            motherp->updateNeighbourBoundary(sigmaneigh, -1);
                            (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(), -1);
                        }
                    } // end neighbour loop
                }
            }
        }
    int expected_area = motherp->Area() + daughterp->Area();
    if (pixel_count != expected_area) {
        throw runtime_error("Cell area is " + to_string(expected_area) + " but only "
                            + to_string(pixel_count) + " pixels were found inside bounding box");
    }

    delete[] (celldir);

    return daughterp->sigma;
}

int CellularPotts::PlaceCellsRandomly(int n, int cellsize) {
    bool overlap = false;
    int radsq = (int) (((double) cellsize) / 3.14);
    int count = 0, check = 0;

    while (count < n && check < 100) {
        overlap = false;
        int x0 = radsq + RandomNumber(sizex - 2 * radsq); //should avoid putting them across boundaries, radsq
        //(radius square is overdoing it,
        //but it's ok given the bugs that happen
        //when you put cells across boundaries :S )
        int y0 = radsq + RandomNumber(sizey - 2 * radsq);
        // check overlap
        for (int x = x0 - cellsize; x < x0 + cellsize; x++) {
            for (int y = y0 - cellsize; y < y0 + cellsize; y++) {
                if ((x - x0) * (x - x0) + (y - y0) * (y - y0) < radsq && x >= 1 && x < sizex - 1 && y >= 1 &&
                    y < sizey - 1)  // circle
//              if( abs(x-x0)<sqrt(radsq) && abs(y-y0)<sqrt(radsq)+5 && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) //rectangle
//             if( ((x-x0)*(x-x0)/((double)radsq-20) +(y-y0)*(y-y0)/(100+(double)radsq))<1. && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)  // ellipse

                {
                    if (sigma[x][y]) {
                        overlap = true;
                        check++;
                        //cerr << "Overlap. count: "<<n<<" check: "<<check<< endl;
                        break;
                    }
                }

            }
            if (overlap)
                break;
        }

        if (!overlap) {
            check = 0;
            count++;
            //cerr << "No overlap. count: "<<n<<" check: "<<check<< endl;
            for (int x = x0 - cellsize; x < x0 + cellsize; x++) {
                for (int y = y0 - cellsize; y < y0 + cellsize; y++) {
//                   if( abs(x-x0)<sqrt(radsq) && abs(y-y0)<sqrt(radsq)+5 && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)   //rectangle
                    if ((x - x0) * (x - x0) + (y - y0) * (y - y0) < radsq && x >= 1 && x < sizex - 1 && y >= 1 &&
                        y < sizey - 1) //circle
//                  if( ((x-x0)*(x-x0)/((double)radsq-20) +(y-y0)*(y-y0)/(100+(double)radsq))<1.  && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) //ellipses
                        sigma[x][y] = count;

                }

            }
        }
    }

    cout << "Placed " << count << " cells out of " << n << " requested." << endl;
    //exit(1);
    return count;
}

// Places cells at regular distance from one another:
// For square cells of size s, the spatial occupation is sqrt(s).
// Since we have to place n of them, we use a square of space of size sqrt(n)*(sqrt(s) + a_little_bit),
// centered at the center of grid, which means that the upper left corner of the first cell is at
// x=(sizex-sqrt(n)*(sqrt(s) + a_little_bit))/2
int CellularPotts::PlaceCellsOrderly(int n_cells, int size_cells) {
    int count = 0;
    int a_little_bit = 2;
    int beginx, endx, beginy, endy;
    int step, avrg_area;
    int smaller_dimension = (par.sizex < par.sizey) ? par.sizex : par.sizey;
    int sqrt_n_cells = 1 + sqrt(n_cells);
    //to avoid having 49 cells when you want 50, I'm rounding sqrt(n_cells) to +1
    if ((sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) > smaller_dimension) {
        std::cerr << "PlaceCellsOrderly(): Error. Too many cells or too large size?" << '\n';
        exit(1);
    }

    step = (sqrt(size_cells) + a_little_bit);

    avrg_area = 0;

    // each x,y point denotes the upper left corner of a cell
    // with i,j we run through the cell
    // for different initial conditions we will have oders of putting the cells
    //this is quite extendable
    vector<int> v_order_x;
    vector<int> v_order_y;

    if (par.cell_placement) {
        std::cerr << "yay, placing cells where I want!" << endl;
        switch (par.cell_placement) {
            case 1:
                beginx = a_little_bit;
                endx = a_little_bit + sqrt_n_cells * (sqrt(size_cells) + a_little_bit);
                beginy = (par.sizey - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                endy = (par.sizey + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                break;
            case 2:
                beginx = (par.sizex - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                endx = (par.sizex + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                beginy = a_little_bit;
                endy = a_little_bit + sqrt_n_cells * (sqrt(size_cells) + a_little_bit);
                break;
            case 3:
                beginx = par.sizex - sqrt_n_cells * (sqrt(size_cells) + a_little_bit) - a_little_bit;
                endx = par.sizex - a_little_bit;
                beginy = (par.sizey - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                endy = (par.sizey + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                break;
            case 4:
                beginx = (par.sizex - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                endx = (par.sizex + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
                beginy = par.sizey - sqrt_n_cells * (sqrt(size_cells) + a_little_bit) - a_little_bit;
                endy = par.sizey - a_little_bit;
                break;
            default:
                std::cerr << "PlaceCellsOrderly(): Error. Got an unusable value for par.cell_placement" << '\n';
                exit(1);
        }
        for (int x = beginx; x < endx; x += step) v_order_x.push_back(x);
        for (int y = beginy; y < endy; y += step) v_order_y.push_back(y);
    } else {
        beginx = (par.sizex - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
        endx = (par.sizex + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
        beginy = (par.sizey - sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;
        endy = (par.sizey + sqrt_n_cells * (sqrt(size_cells) + a_little_bit)) / 2;


        switch (par.init_cell_config) {
            case 0: {
                for (int x = beginx; x < endx; x += step) v_order_x.push_back(x);
                for (int y = beginy; y < endy; y += step) v_order_y.push_back(y);
            }
                break;
            case 1: {
                for (int x = endx; x > beginx; x -= step) v_order_x.push_back(x);
                for (int y = beginy; y < endy; y += step) v_order_y.push_back(y);
            }
                break;
            case 2: {
                for (int x = beginx; x < endx; x += step) v_order_x.push_back(x);
                for (int y = endy; y > beginy; y -= step) v_order_y.push_back(y);
            }
                break;
            case 3: {
                for (int x = endx; x > beginx; x -= step) v_order_x.push_back(x);
                for (int y = endy; y > beginy; y -= step) v_order_y.push_back(y);
            }
                break;
            default:
                std::cerr << "PlaceCellsOrderly(): Error. Got an unusable value for par.init_config" << '\n';
                exit(1);
        }
    }

    for (auto x: v_order_x) {
        for (auto y: v_order_y) {
            // for(int x = beginx ; x < endx ; x += step ){
            //   for(int y = beginy ; y < endy ; y += step ){
            std::cout << "Cell will be placed at: " << x << "," << y << '\n';

            count++;
            int this_area = 0;
            for (int i = 0; i < sqrt(size_cells); i++) {
                for (int j = 0; j < sqrt(size_cells); j++) {
                    if (sigma[x + i][y + j]) {
                        std::cerr << "Sigma = " << sigma[x + i][y + j] << '\n';
                        std::cerr << "Grid point " << x + i << "," << y + j << " is already occupied" << '\n';
                        exit(1);
                    }
                    sigma[x + i][y + j] = count;
                    this_area++;
                    avrg_area++;
                    if (this_area == size_cells) break;
                }
                if (this_area == size_cells) break;
            }
            if (count == n_cells) break;
        }
        if (count == n_cells) break;
    }


    cerr << "Placed " << count << " cells out of " << n_cells << " requested; avrg area = "
         << avrg_area / (double) count
         << endl;
    //exit(1);
    return count;
}

void CellularPotts::killCell(int c_sigma) {
    Cell &c = (*cell)[c_sigma];
    c.SetTargetArea(0);
    c.Apoptose(); //set alive to false
    RemoveCell(&c, par.min_area_for_life, int(c.meanx), int(c.meany));
}

//WARNING: Cells placed with this function will be much bigger than cell_size: algoritm doesn't work properly
int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y) {

    // make initial cells using Eden Growth

    int **new_sigma = (int **) malloc(sizex * sizeof(int *));
    if (new_sigma == nullptr)
        MemoryWarning();

    new_sigma[0] = (int *) malloc(sizex * sizey * sizeof(int));
    if (new_sigma[0] == nullptr)
        MemoryWarning();

    for (int i = 1; i < sizex; i++)
        new_sigma[i] = new_sigma[i - 1] + sizey;

    /* Clear CA plane */
    {
        for (int i = 0; i < sizex * sizey; i++)
            new_sigma[0][i] = 0;
    }


    // scatter initial points, or place a cell in the middle
    // if only one cell is desired
    int cellnum = cell->size() - 1;

    if (n_cells > 1) {


        {
            for (int i = 0; i < n_cells; i++) {

                sigma[RandomNumber(sx) + offset_x][RandomNumber(sy) + offset_y] = ++cellnum;

            }
        }
    } else {
        sigma[sx][sy] = ++cellnum;

    }

    // Do Eden growth for a number of time steps
    {
        for (int i = 0; i < cell_size; i++) {
            for (int x = 1; x < sizex - 1; x++)
                for (int y = 1; y < sizey - 1; y++) {

                    if (sigma[x][y] == 0) {
                        // take a random neighbour
                        int xyp = (int) (8 * RANDOM() + 1);
                        int xp = nx[xyp] + x;
                        int yp = ny[xyp] + y;
                        int kp;
                        //  NB removing this border test yields interesting effects :-)
                        // You get a ragged border, which you may like!
                        if ((kp = sigma[xp][yp]) != -1)
                            if (kp > (cellnum - n_cells))
                                new_sigma[x][y] = kp;
                            else
                                new_sigma[x][y] = 0;
                        else
                            new_sigma[x][y] = 0;

                    } else {
                        new_sigma[x][y] = sigma[x][y];
                    }
                }

            // copy getSigma to new_sigma, but do not touch the border!
            {
                for (int x = 1; x < sizex - 1; x++) {
                    for (int y = 1; y < sizey - 1; y++) {
                        sigma[x][y] = new_sigma[x][y];
                    }
                }
            }
        }
    }
    free(new_sigma[0]);
    free(new_sigma);

    return cellnum;
}

int CellularPotts::FancyloopX(int loopdepth, int meanx, int meany, int thissig, bool above) {
    int removed = 0;
    bool loop = true;
    int py, px;
    int sigmaneigh;
    py = meany - ((int) above * 2 - 1) * loopdepth;
    if (py <= 0 || py >= sizey - 1) {
        if (par.periodic_boundaries) {
            if (py <= 0)
                py = sizey - 2 + py;
            else if (py >= sizey - 1)
                py = py - sizey + 2;
        } else {
            loop = false;
        }
    }

    if (loop)
        for (int x = meanx - loopdepth; x <= meanx + loopdepth; x++) {
            px = x;
            if (x <= 0 || x >= sizex - 1) {
                if (par.periodic_boundaries) {
                    if (x <= 0)
                        px = sizex - 2 + x;
                    else if (x >= sizex - 1)
                        px = x - sizex + 2;
                } else {
                    continue;
                }
            }
//         cerr<<"removing?"<<endl;
            if (sigma[px][py] == thissig) {
//           cerr<<"going on with removing X"<<endl;
                sigma[px][py] = 0;
                removed++;
                for (int k = 1; k <= n_nb; k++) {
                    int neix = px + nx[k];
                    int neiy = py + ny[k];
                    if (neix <= 0 || neix >= sizex - 1 || neiy <= 0 || neiy >= sizey - 1) {
                        if (par.periodic_boundaries) {
                            if (neix <= 0) neix = sizex - 2 + neix;
                            if (neix >= sizex - 1) neix = neix - sizex + 2;
                            if (neiy <= 0) neiy = sizey - 2 + neiy;
                            if (neiy >= sizey - 1) neiy = neiy - sizey + 2;
                        } else {
                            continue;
                        }
                    }
                    sigmaneigh = sigma[neix][neiy];
                    if (sigmaneigh == 0) {
                        try {
                            edgeSetVector.erase({px, py, neix, neiy});
                            edgeSetVector.erase({neix, neiy, px, py});
                        } catch (const std::out_of_range &oor) {}
                    }
                }
            }
        }

    return removed;
}

int CellularPotts::FancyloopY(int loopdepth, int meanx, int meany, int thissig, bool left) {
    int removed = 0;
    bool loop = true;
    int py, px;
    int sigmaneigh;

    px = meanx - ((int) left * 2 - 1) * loopdepth;
    if (px <= 0 || px >= sizex - 1) {
        if (par.periodic_boundaries) {
            if (px <= 0)
                px = sizex - 2 + px;
            else if (px >= sizex - 1)
                px = px - sizex + 2;
        } else {
            loop = false;
        }
    }

    if (loop)
        for (int y = meany - loopdepth + 1; y <= meany + loopdepth - 1; y++) {
            py = y;
            if (y <= 0 || y >= sizey - 1) {
                if (par.periodic_boundaries) {
                    if (y <= 0)
                        py = sizey - 2 + y;
                    else if (y >= sizey - 1)
                        py = y - sizey + 2;
                } else {
                    continue;
                }
            }

            if (sigma[px][py] == thissig) {
                sigma[px][py] = 0;
                removed++;
                for (int k = 1; k <= n_nb; k++) {
                    int neix = px + nx[k];
                    int neiy = py + ny[k];
                    if (neix <= 0 || neix >= sizex - 1 || neiy <= 0 || neiy >= sizey - 1) {
                        if (par.periodic_boundaries) {
                            if (neix <= 0) neix = sizex - 2 + neix;
                            if (neix >= sizex - 1) neix = neix - sizex + 2;
                            if (neiy <= 0) neiy = sizey - 2 + neiy;
                            if (neiy >= sizey - 1) neiy = neiy - sizey + 2;
                        } else {
                            continue;
                        }
                    }
                    sigmaneigh = sigma[neix][neiy];
                    if (sigmaneigh == 0) {
                        try {
                            edgeSetVector.erase({px, py, neix, neiy});
                            edgeSetVector.erase({neix, neiy, px, py});
                        } catch (const std::out_of_range &oor) {}
                    }
                }
            }
        }

    return removed;
}


//Hopefully a little less expensive function than the one above
void CellularPotts::RemoveCell(Cell *thiscell, int min_area, int meanx, int meany) {
    int sigmaneigh, thissig, thisarea;
    int countpix = 0; //to check if we removed all pixels
    bool loop = true;
    int loopdepth = 1;
    thissig = thiscell->Sigma();
    thisarea = thiscell->Area();

//   cerr<<"meanx: "<<meanx<<" meany: "<<meany<<endl;
    if (!thisarea) {
        cerr << "CA-RemoveCell warning: attempting to remove dead cell. " << endl;
    }

    if (sigma[meanx][meany] == thissig) {
        sigma[meanx][meany] = 0;
        countpix++;
//     cerr<<"Removed one, countpix ="<<countpix<<endl;
    }

    // TODO: Profile this version against looping over BoundingBox
    while (true) {
        //top row
        //Here we go through progressively larger squares (only its boundary)
        //first through the top row, then the bottom row, then left and right sides
        //which are smaller.
        countpix += FancyloopX(loopdepth, meanx, meany, thissig, true);
        if (countpix == thisarea) break;
        countpix += FancyloopX(loopdepth, meanx, meany, thissig, false);
        if (countpix == thisarea) break;
        countpix += FancyloopY(loopdepth, meanx, meany, thissig, true);
        if (countpix == thisarea) break;
        countpix += FancyloopY(loopdepth, meanx, meany, thissig, false);
        if (countpix == thisarea) break;

//     cerr<<"Area: "<<thisarea<<" loopdepth: "<<loopdepth<<" countpix: "<<countpix<<endl;
//     if(loopdepth>10) exit(1);

        loopdepth++;

    }

    //take care that area is set to 0
    thiscell->DecrementAreaBy(thisarea);
    //removed the tracking of the moments. Check if this gives major problems


    //deal with the cells that have this cell as a neighbour
    for (auto &neigh: thiscell->neighbours) {
        int signeigh = neigh.first;
        int blength = neigh.second.first;

        if (signeigh) {
            (*cell)[signeigh].setNeighbour(thissig, 0, 0);
            (*cell)[signeigh].updateNeighbourBoundary(0, blength);
        }

    }

    thiscell->clearNeighbours();

//   exit(1);

}

// Predicate returns true when connectivity is locally preserved
// if the value of the central site would be changed
// note that at this point we already know that focal point is neighbouring a point with different getSigma
bool CellularPotts::ConnectivityPreservedP(int x, int y) {

    // Use local nx and ny in a cyclic order (starts at upper left corner)
    // zeroth site is ignored and first site is repeated, for easier looping - that's why the array is 10 long
    // even though there are only 8 neighbours
    // (see below: int s_next_nb=getSigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]];)
    const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1};
    const int cyc_ny[10] = {0, -1, -1, -1, 0, 1, 1, 1, 0, -1};

    int sxy = sigma[x][y]; // the central site
    if (sxy == 0) return true;

    int n_borders = 0; // to count the amount of sites in state sxy bordering a site !=sxy

    static int stack[8]; // stack to count number of different surrounding cells
    int stackp = -1;
    bool one_of_neighbors_medium = false;

    for (int i = 1; i <= 8; i++) {

        int s_nb = sigma[x + cyc_nx[i]][y + cyc_ny[i]]; //this is getSigma of neighbour
        int s_next_nb = sigma[x + cyc_nx[i + 1]][y + cyc_ny[i + 1]]; //this is getSigma of next neighbour on the list
        //if at least one of them == focal getSigma, and they are different from each other
        // i.e. if one of them is same as sxy and other not (but whihc one does not matter)
        if ((s_nb == sxy || s_next_nb == sxy) && (s_nb != s_next_nb)) {

            // check whether s_nb is adjacent to non-identical site,
            n_borders++; // count it
        }
        int j;
        bool on_stack_p = false;

        // we need the next heuristic to prevent stalling at
        // cell-cell borders
        // do not enforce constraint at two cell interface(no medium)
        if (s_nb) {
            for (j = stackp; j >= 0; j--) {
                if (s_nb == stack[j]) {
                    on_stack_p = true;
                    break;
                }
            }
            if (!on_stack_p) {
                if (stackp > 6) {
                    cerr << "Stack overflow, stackp=" << stackp << "\n";
                }
                stack[++stackp] = s_nb;
            }
        } else {
            one_of_neighbors_medium = true;
        }
    }

    // number of different neighbours is stackp+1;
    if (n_borders > 2 && ((stackp + 1) > 2 || one_of_neighbors_medium)) {
        return false;
    } else
        return true;

}


