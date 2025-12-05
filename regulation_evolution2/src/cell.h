/*

Copyright 1996-2006 Roeland Merks

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
#ifndef _CELL_HH_
#define _CELL_HH_

#include "parameter.h"
//#define EMPTY -1
#include <math.h>
#include <map>
#include <utility>
#include "random.h"
#include "genome.h"
#include "boundingbox.h"
#include "parameter.h"

#define MIGRATE 1
#define DIVIDE 2

class Dish;

class Cell {
    friend class Dish;
    friend class CellularPotts;

public:


    /*! \brief Constructor to insert a cell into Dish "who"

    Used to add a new Cell to the dish: new Cell(dish,
    celtype).
    */
    explicit Cell(const Dish &who, int settau = 1, int setrecycledsigma = -1) {
        owner = &who;

        alive = true;
        group = par.randgroups ? RANDOM() < 0.5 : 0;

        meanx = 0.;
        meany = 0.;
        chemmu = 0.0;
        dividecounter = 0;
        grad_conc = 0;

        colour_of_birth = 1;
        date_of_birth = 0;
        times_divided = 0;
        mother = 0;
        daughter = 0;

        // add new elements to each of the dimensions of "J"


        // maxsigma keeps track of the last cell identity number given out to a cell
        if (setrecycledsigma == -1) {
            // amount gives the total number of Cell instantiations (including copies)
            amount++;

            sigma = maxsigma++;
            //cout<<"check: not recycling, new sigma is "<<getSigma<<endl;
        } else {
            sigma = setrecycledsigma;
            //cout<<"check: YES recycling, new sigma is "<<getSigma<<endl;
        }

        ancestor = sigma;

        //if (!J) {
        //  ReadStaticJTable(par.Jtable);
        //}

        // This should not be here
        //if(!vJ){
        //  ReadKeyLockFromFile(par.keylock_list_filename)
        //}

        //create the genome: use parameters for size. This is a randomly generated genome.
        //genome.InitGenome(par.innodes, par.regnodes, par.outnodes);
        //let's do this separately.
        gextiming = 0;
        tau = settau;
        area = 0;
        target_area = 0;
        length = 0;
        sum_x = 0;
        sum_y = 0;
        sum_xx = 0;
        sum_yy = 0;
        sum_xy = 0;

        v[0] = 0.;
        v[1] = 0.;
        n_copies = 0;
        mu = 0.0;
        tvecx = 0.;
        tvecy = 0.;
        prevx = 0.;
        prevy = 0.;

        chemvecx = 0.;
        chemvecy = 0.;

        food = par.foodstart;
        last_meal = 0;

        time_since_birth = 0;

        persdur = 0;
        perstime = 0;
        if (par.n_chem) {
            chem = new double[par.n_chem];
        }
    }

    ~Cell();


    //! Default copy constructor.
    // Copies a cell into a new one
    Cell(const Cell &src) {
        //cout<<"Tomato"<<endl;
        // make an exact copy (for internal use)
        sigma = src.sigma;
        amount++;
        area = src.area;
        target_area = src.target_area;
        length = src.length;
        target_area = src.target_area;

        ancestor = src.ancestor;
        mother = src.mother;
        daughter = src.daughter;
        times_divided = src.times_divided;
        date_of_birth = src.date_of_birth;
        //cerr<<"this?: "<<date_of_birth<<endl;
        colour_of_birth = src.colour_of_birth;
        group = src.group;
        dividecounter = src.dividecounter;
        genome = src.genome;
        gextiming = src.gextiming;

        tau = src.tau;
        alive = src.alive;
        v[0] = src.v[0];
        v[1] = src.v[1];
        n_copies = src.n_copies;
        sum_x = src.sum_x;
        sum_y = src.sum_y;
        sum_xx = src.sum_xx;
        sum_yy = src.sum_yy;
        sum_xy = src.sum_xy;

        meanx = src.meanx;
        meany = src.meany;
        tvecx = src.tvecx;
        tvecy = src.tvecy;
        prevx = src.prevx;
        prevy = src.prevy;
        persdur = src.persdur;
        perstime = src.perstime;
        mu = src.mu;

        chemmu = src.chemmu;
        chemvecx = src.chemvecx;
        chemvecy = src.chemvecy;
        grad_conc = src.grad_conc;

        food = src.food;
        last_meal = src.last_meal;

        owner = src.owner;

        neighbours = src.neighbours;
        jkey_dec = src.jkey_dec;
        jlock_dec = src.jlock_dec;

        if (par.n_chem) {
            chem = new double[par.n_chem];
            for (int ch = 0; ch < par.n_chem; ch++)
                chem[ch] = src.chem[ch];
        }

        time_since_birth = src.time_since_birth;
    }

    /*! \brief Add a new cell to the dish.

       Call it as: new Cell(parent, true); mother will be modified for
       ancestry administration!

       \param settau
       Cell type of daughter cell.
    */
    void CellBirth(Cell &mother);

    /*! \brief Assignment operator.

    Called if one cell is assigned to another. Remember to change both
    assignment operator and copy constructor when adding new attributes
    to Cell.
    */
    inline Cell &operator=(const Cell &src) {
        if (this == &src)
            return *this;
        //cout<<"Potato"<<endl;
        group = src.group;
        alive = src.alive;
        sigma = src.sigma;
        area = src.area;
        tau = src.tau;
        target_area = src.target_area;
        date_of_birth = src.date_of_birth;
        times_divided = src.times_divided;
        genome = src.genome;
        gextiming = src.gextiming;
        dividecounter = src.dividecounter;

        v[0] = src.v[0];
        v[1] = src.v[1];
        n_copies = src.n_copies;

        meanx = src.meanx;
        meany = src.meany;
        tvecx = src.tvecx;
        tvecy = src.tvecy;
        prevx = src.prevx;
        prevy = src.prevy;

        persdur = src.persdur;
        perstime = src.perstime;
        mu = src.mu;

        chemmu = src.chemmu;
        chemvecx = src.chemvecx;
        chemvecy = src.chemvecy;

        food = src.food;
        last_meal = src.last_meal;

        grad_conc = src.grad_conc;
        sum_x = src.sum_x;
        sum_y = src.sum_y;
        sum_xx = src.sum_xx;
        sum_yy = src.sum_yy;
        sum_xy = src.sum_xy;

        length = src.length;
        target_area = src.target_area;
        amount++;
        owner = src.owner;

        neighbours = src.neighbours;

        ancestor = src.ancestor;

        chem = new double[par.n_chem];
        for (int ch = 0; ch < par.n_chem; ch++)
            chem[ch] = src.chem[ch];

        time_since_birth = src.time_since_birth;

        return *this;

    }

    /*! \brief Returns false if Cell has apoptosed (vanished). */
    inline bool AliveP() const {
        return alive;
    }

    //! Returns the cell colour.
    inline int getGroup() const {
        return group;
    };

//this function maps migration vector angle to a colour in radial_colour array (see misc.cpp)
    inline int AngleColour() const {
        double ang = atan(tvecy / tvecx);
        if (tvecx < 0.000) ang += M_PI;
        else if (tvecy < 0.000) ang += 2 * M_PI;
        ang /= 2. * M_PI;

        return (int) (ang * 254) + 2;
    };

    //sets properties of cell
    // for now only growth   ---- THIS STUFF CREATES PROBLEMS!

    inline void InitMeanX(double xpos) {
        meanx = xpos;
    }

    inline void InitMeanY(double ypos) {
        meany = ypos;
    }

    //Return values related to cell position and direction of migration
    inline double getXpos() const {
        return meanx;
    }

    inline double getYpos() const {
        return meany;
    }

    inline double getXvec() const {
        return tvecx;
    }

    inline double getYvec() const {
        return tvecy;
    }

    inline double getChemXvec() const {
        return chemvecx;
    }

    inline double getChemYvec() const {
        return chemvecy;
    }

    inline double getChemMu() const {
        //cout<<"numu: "<<mu<<endl;
        return chemmu;
    }

    inline double getMu() const {
        //cout<<"numu: "<<mu<<endl;
        return mu;
    }

    inline void setMu(double initmu) {
        mu = initmu;
        // cout<<"initmu: "<<mu<<endl;
    }

    inline void setChemMu(double initweight) {
        chemmu = initweight;
        // cout<<"initmu: "<<mu<<endl;
    }

    inline void setPersDur(int dur) {
        persdur = dur;
    }

    inline void setPersTime(int time) {
        perstime = time;
    }

    //initialise variables for target vector position
    inline void startTarVec() {
        prevx = meanx;
        prevy = meany;
        double pol = RANDOM() * 2. * M_PI; //random angle

        //pol=0.;

        tvecx = cos(pol);  // try swapping these around
        tvecy = sin(pol);
    }

    inline void startChemVec() { //to make sure hamiltonian has something to work with.

        double pol = RANDOM() * 2. * M_PI; //random angle

        //pol=0.;

        chemvecx = cos(pol);  // try swapping these around
        chemvecy = sin(pol);
    }

    inline void setChemVec(double xx, double yy) { //to make sure hamiltonian has something to work with.
        chemvecx = xx;  // try swapping these around
        chemvecy = yy;
    }

    //update the target vector with the actual direction of motion
    inline void updateTarVec() {

        if ((meanx - prevx) * (meanx - prevx) > 0.0000001 || (meany - prevy) * (meany - prevy) > 0.0000001) {
            tvecx = meanx - prevx;
            tvecy = meany - prevy;
            double hyphyp = hypot(tvecx, tvecy);
            tvecx /= hyphyp;
            tvecy /= hyphyp;
            prevx = meanx;
            prevy = meany;
        } else {
            //just keep moving in the same direction then
            prevx = meanx;
            prevy = meany;
        }
    }

    inline void updatePersTime() {
        perstime++;
        if (perstime == persdur) {
            updateTarVec();
            //cerr<<"persdur="<<persdur<<endl;
            //cerr<<"Cell "<<getSigma<<" geupdate, "<<tvecx<<" "<<tvecy<<endl;
            perstime = 0;
        }
    }

    double MAXmu = 30;

    inline void MutateMu() {
        mu += (RANDOM() - 0.5) / 1.;
        if (mu < 0) mu = -mu;
        else if (mu > MAXmu) mu = MAXmu - mu;
    }

    //! Set cell type of this Cell.
    inline void setTau(int settau) {
        tau = settau;
    }

    //! Get cell type of this Cell.
    inline int getTau() const {
        return tau;
    }

    //! Set group of this cell to new_group, irrespective of type.
    inline void setGroup(const int new_group) {
        group = new_group;
    }

    //! Return Cell's actual area.
    inline int Area() const {
        return area;
    }

    //! Return Cell's target area.
    inline int TargetArea() const {
        return target_area;
    }


    //! Debugging function used to print the cell's current inertia tensor (as used for calculations of the length )
    inline void PrintInertia() const {

        double ixx = (double) sum_xx - (double) sum_x * sum_x / (double) area;
        double iyy = (double) sum_yy - (double) sum_y * sum_y / (double) area;
        double ixy = (double) sum_xy - (double) sum_x * sum_y / (double) area;

        cout << "ixx = " << ixx << "\n";
        cout << "iyy = " << iyy << "\n";
        cout << "ixy = " << ixy << "\n";

    }

    // return the current length
    inline double Length() const {
        return length;
    }

    /*! \brief Returns the maximum cell identity number in the Dish.
      This would normally be the number of cells in the Dish, although
     the number includes apoptosed cells.
    */
    static inline int MaxSigma() {
        return maxsigma;
    }

    //! Returns the cell's cell identity number.
    inline int Sigma() const {
        return sigma;
    }

    // THIS REALLY DOES NOT WORK, because getSigma is protected.
//   // set new getSigma at birth
//   inline int SetSigmaIfRecycled(const int recycled_sigma) const {
//     return getSigma=recycled_sigma;
//   }

    //! Sets the target area of the cell.
    inline int SetTargetArea(const int new_area) {
        return target_area = new_area;
    }

    //! Sends the current cell into apoptosis
    inline void Apoptose() {
        alive = false;
    }

    //! Decrement the cell's target area by one unit.
    inline int IncrementTargetArea() {
        return ++target_area;
    }

    //! Increment the cell's target area by one unit.
    inline int DecrementTargetArea() {
        return --target_area;
    }

    //! This is the oldest ancestor of this cell since last time we saved the ancestry
    inline int getAncestor() const { return ancestor; }

    //! Resets the ancestor to the cell itself. Called whenever the ancestry information is saved.
    inline void resetAncestor() {
        ancestor = Sigma();
    }

    //! Cell lineage tracking: get the cell's parent
    inline int Mother() const { return mother; }

    //! Cell lineage tracking: get the cell's daughter
    inline int Daughter() const { return daughter; }

    //! Returns a counter keeping track of the number of divisions
    inline int TimesDivided() const { return times_divided; }

    inline void AddTimesDivided() { times_divided++; }

    inline void ResetTimesDivided() {
        times_divided = 0;
        dividecounter = 0;
    }

    //! Returns Monte Carlo Step (MCS) when this cell originated.
    inline int DateOfBirth() const { return date_of_birth; }

    //! Returns the cell type at the time of birth.
    inline int ColourOfBirth() const { return colour_of_birth; }

    //! Returns the bond energy J between this cell and cell c2.
    inline int GetJ(const Cell &c2) const {
        return J[sigma][c2.sigma];
    }


    //! Sets bond energy J between cell type t1 and t2 to val
    inline static int SetJ(int t1, int t2, int val) {
        return J[t2][t1] = J[t1][t2] = val;
    }


    //deal with the genome:
    inline void ClearGenomeState() {
        genome.ResetGenomeState();
    }

    inline void ReadGenomeFromFile(char *filename) {
        genome.ReadFromFile(filename);
    }

    inline void WriteGenomeToFile(char *filename) {
        genome.WriteToFile(filename);
    }

    inline void CreateRandomGenome(int in, int reg, int out, vector<double> const &inscale) {
        genome.InitGenome(in, reg, out, inscale);
    }

    inline void UpdateGenes(const array<double, 2> &input, bool sync) {
        genome.UpdateGeneExpression(input, sync);
    }

    inline void FinishGeneUpdate() {
        genome.FinishUpdate();
    }

    inline void MutateGenome(double mu_, double mustd) {
        genome.MutateGenome(mu_, mustd);
    }

    inline void setGTiming(int timing) {
        gextiming = timing;
    }

    inline int Gextiming() const {
        return gextiming;
    }

    //This is for keeping track of who is neigh, how much contact and for how long
    std::map<int, pair<int, int>> neighbours; //stores neigh as {neighbouring cells(ID) <amount of contact,duration>
    //cell neighbours
    void setNeighbour(int neighbour, int boundarylength, int contactduration);

    void clearNeighbours();

    int returnBoundaryLength(int cell);

    int returnDuration(int cell);

    int SetNeighbourDurationFromMother(int cell, int motherduration);

    int updateNeighbourBoundary(int cell, int boundarymodification);

    int updateNeighbourDuration(int cell, int durationmodification);

    // Returns a bounding box around the cell where CPM->getSigma(x, y) is likely to be cell->sigma
    // Instead of dynamically calculating this, we could keep it is an attr and update it on ConvertSpin
    // Is that faster?
    BoundingBox getBoundingBox();

    void updateJDecs();

private:
    // used internally by class CellularPotts
    inline void CleanMoments() {
        sum_x = sum_y = sum_xx = sum_xy = sum_yy = area = target_area = 0;
    }

    // used internally by class CellularPotts
    // notice that it should be preceded always by IncrementArea()
    //
    // so far: problem is I am not changing sumx, only meanx


    inline double AddSiteToMoments(int x, int y, double new_l = -1.) {

        //calculate higher moments ONLY if not wrapped boundaries
        length = 0.; //set length to zero becaues we don't need it anyway

        // Add a site to the raw moments, then update and return the
        // length of the cell
        //cout<<"area now: "<<area<<"x: "<<x<<"y: "<<y<<endl;
        // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
        // Eventually this function may be used to carry
        // out all necessary adminstration at once
        if (par.periodic_boundaries) {
            if (area > 0) {
                //also notice that this is very related to meanx...

                //double curr_avrg_x = sum_x/((double)(area-1)); //area-1 because we have just incremented area without counting this
                //double curr_avrg_y = sum_y/((double)(area-1));

                //if x is closer to running average when wrapped, we wrap it
                if ((x - meanx) > 0 && (x - meanx) > (meanx - (x - (par.sizex - 2)))) {
                    //if( abs(meanx - x) > abs(meanx - (x - (par.sizex-2))) ){
                    //print passing border
                    //cerr<<"Astm: meanx ="<<meanx<<", pixel on the right: "<<x<<", we wrap it to "<<x-(par.sizex-2)<<endl;
                    x -= (par.sizex -
                          2); //par.sizex -2 because to compare floats with int we assume the pixel begins (on the left, with integer)
                    //cerr<<"passb1"<<endl;
                    //cerr<<"so sum_x will be: "<<sum_x+x<<endl;
                } else if ((meanx - x) > 0 && (meanx - x) > (x + (par.sizex - 2) - meanx)) {
                    //else if( abs(meanx - x) > abs(meanx - (x + (par.sizex-2))) ) {
                    //print passing border
                    //cerr<<"Astm: meanx ="<<meanx<<", pixel on the left: "<<x<<", we wrap it to "<<x+(par.sizex-2)<<endl;
                    x += (par.sizex - 2);
                    //cerr<<"so sum_x will be: "<<sum_x+x<<endl;
                    //cerr<<"passb2"<<endl;
                }
                //same for y
                if ((y - meany) > 0 && (y - meany) > (meany - (y - (par.sizey - 2)))) {
                    y -= (par.sizey - 2);
                    //cerr<<"passb3 meany: "<<meany<<endl;
                } else if ((meany - y > 0) && meany - y > (y + (par.sizey - 2) - meany)) {
                    y += (par.sizey - 2);
                    //cerr<<"passb4"<<endl;
                }


                sum_x += x;
                //cerr<<"Astm sum_x= "<<sum_x<<endl;
                sum_y += y;

                //cerr<<"Before add: meanx: "<<meanx<<", meany: "<<meany<<endl;

                //now if meanx or meany are outside borders we wrap them
                // and shift prevx or y accordingly
                meanx = sum_x / ((double) (area));
                meany = sum_y / ((double) (area));



                //if x is closer to running average when wrapped, we wrap it


                //notice that there is equal sign with if( meanx >= sizex-1 )
                // this means that the right edge of the CA does not exist
                // (it is the first point wrapped)
                if (meanx >= par.sizex - 1) {
                    //wrap meanx
                    //cerr<<"meanx>par.sizex-1 we wrap meanx, sum_x and prevx"<<endl;
                    meanx -= (par.sizex - 2);
                    //change sumx and prevx
                    sum_x -= area * (par.sizex -
                                     2); // sumx has to be shifted, as a whole (each pixel by sizex-2, so in total by area*(sizex-2))
                    prevx -= (par.sizex - 2);
                    //cerr<<"passb5"<<endl;
                } else if (meanx < 1) {
//           cerr<<"meanx<1 we wrap it meanx, sum_x and prevx"<<endl;
                    meanx += (par.sizex - 2);
                    sum_x += area * (par.sizex - 2);
                    prevx += (par.sizex - 2);
                    //cerr<<"passb6"<<endl;
                }
                //same for y
                if (meany >= par.sizey - 1) {
                    meany -= (par.sizey - 2);
                    sum_y -= area * (par.sizey - 2);
                    prevy -= (par.sizey - 2);
                    //cerr<<"passb7"<<endl;
                } else if (meany < 1) {
                    meany += (par.sizey - 2);
                    sum_y += area * (par.sizey - 2);
                    prevy += (par.sizey - 2);
                    //cerr<<"passb8 new meany: "<<meany<<endl;
                }

            } else {
                cerr << "AddSiteToMoments(): Error. How can area be zero after incrementing it?" << endl;
                exit(1);
            }
            //cerr<<"After add: meanx: "<<meanx<<", meany: "<<meany<<endl;
        }
            //else if not periodic_boundaries
        else {
            //update alll moments
            sum_x += x;
            sum_y += y;
            sum_xx += x * x;
            sum_yy += y * y;
            sum_xy += x * y;

            // update length (see appendix. A, Zajac.jtb03), if length is not given
            // NB. 24 NOV 2004. Found mistake in Zajac's paper. See remarks in
            // method "Length(..)".
            if (new_l < 0.) {
                length = Length(sum_x, sum_y, sum_xx, sum_yy, sum_xy, area);
            } else {
                length = new_l;
            }

            //mean position of cell
            if (area > 0) {
                meany = double(sum_y) / double(area);
                meanx = double(sum_x) / double(area);
            }
        }
        return length;
    }

    // used internally by class CellularPotts
    inline double RemoveSiteFromMoments(int x, int y, double new_l = -1.) {

        length = 0.;

        if (par.periodic_boundaries) {
            if (area > 0) {
                //if x is closer to running average when wrapped, we wrap it
                if ((x - meanx) > 0 && (x - meanx) > (meanx - (x - (par.sizex - 2)))) {
                    //if( abs(meanx - x) > abs(meanx - (x - (par.sizex-2))) ){
                    //print passing border
//           cerr<<"Rstm: meanx ="<<meanx<<", pixel on the right: "<<x<<", we wrap it to "<<x-(par.sizex-2)<<endl;
                    x -= (par.sizex -
                          2); //par.sizex -2 because to compare floats with int we assume the pixel begins (on the left, with integer)
                    //cerr<<"so sum_x will be: "<<sum_x-x<<endl;
                    //cerr<<"passb"<<endl;
                } else if ((meanx - x) > 0 && (meanx - x) > (x + (par.sizex - 2) - meanx)) {
                    //else if( abs(meanx - x) > abs(meanx - (x + (par.sizex-2))) ) {
                    //print passing border
//           cerr<<"Rstm: meanx ="<<meanx<<", pixel on the left: "<<x<<", we wrap it to "<<x+(par.sizex-2)<<endl;
                    x += (par.sizex - 2);
                    //cerr<<"so sum_x will be: "<<sum_x-x<<endl;
                    //cerr<<"passb"<<endl;
                }
                //same for y
                if ((y - meany) > 0 && (y - meany) > (meany - (y - (par.sizey - 2)))) {
                    y -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if ((meany - y > 0) && meany - y > (y + (par.sizey - 2) - meany)) {
                    y += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }

                sum_x -= x;
//         cerr<<"Rstm sum_x= "<<sum_x<<endl;
                sum_y -= y;

                //cerr<<"Before rem: meanx: "<<meanx<<", meany: "<<meany<<endl;

                //now if meanx or meany are outside borders we wrap them
                // and shift prevx or y accordingly
                meanx = sum_x / ((double) (area));
                meany = sum_y / ((double) (area));

                //if x is closer to running average when wrapped, we wrap it


                //this is still to correct for -1 or -2
                if (meanx >= par.sizex - 1) {
                    //wrap meanx
//           cerr<<"meanx>par.sizex-1 we wrap meanx, sum_x and prevx"<<endl;
                    meanx -= (par.sizex - 2);
                    //change sumx and prevx
                    sum_x -= area * (par.sizex -
                                     2); // sumx has to be shifted, as a whole (each pixel by sizex-2, so in total by area*(sizex-2))
                    prevx -= (par.sizex - 2);
                    //cerr<<"passb"<<endl;
                } else if (meanx < 1) {
//           cerr<<"meanx<1 we wrap it meanx, sum_x and prevx"<<endl;
                    meanx += (par.sizex - 2);
                    sum_x += area * (par.sizex - 2);
                    prevx += (par.sizex - 2);
                    //cerr<<"passb"<<endl;
                }
                //same for y
                if (meany >= par.sizey - 1) {
                    meany -= (par.sizey - 2);
                    sum_y -= area * (par.sizey - 2);
                    prevy -= (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                } else if (meany < 1) {
                    meany += (par.sizey - 2);
                    sum_y += area * (par.sizey - 2);
                    prevy += (par.sizey - 2);
                    //cerr<<"passb"<<endl;
                }

            }
            //cerr<<"After rem: meanx: "<<meanx<<", meany: "<<meany<<endl;
        }
            //else if not periodic_boundaries
        else {

            // Remove a site from the raw moments, then update and return the
            // length of the cell

            // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
            // Eventually this function may be used to carry
            // out all necessary adminstration at once
            sum_x -= x;
            sum_y -= y;


            sum_xx -= x * x;
            sum_yy -= y * y;
            sum_xy -= x * y;

            // update length (see app. A, Zajac.jtb03), if length is not given
            if (new_l < 0.) {
                length = Length(sum_x, sum_y, sum_xx, sum_yy, sum_xy, area);
            } else {
                length = new_l;
            }

            if (area > 0) {
                meany = double(sum_y) / double(area);
                meanx = double(sum_x) / double(area);
            }
        }
        return length;
    }


    //! \brief Calculates the length based on the given inertia tensor
    //components (used internally)
    static inline double Length(long int s_x, long int s_y, long int s_xx,
                                long int s_yy, long int s_xy, long int n) {

        // inertia tensor (constructed from the raw momenta, see notebook)
        double iyy = (double) s_xx - (double) s_x * s_x / (double) n;
        double ixx = (double) s_yy - (double) s_y * s_y / (double) n;
        double ixy = -(double) s_xy + (double) s_x * s_y / (double) n;

        double rhs1 = (ixx + iyy) / 2., rhs2 = sqrt((ixx - iyy) * (ixx - iyy) + 4 * ixy * ixy) / 2.;

        double lambda_b = rhs1 + rhs2;
        //double lambda_a=rhs1-rhs2;

        // according to Zajac et al. 2003:
        //return 2*sqrt(lambda_b);
        // Grumble, this is not right!!!
        // Must divide by mass!!!!!!

        // see: http://scienceworld.wolfram.com/physics/MomentofInertiaEllipse.html
        //    cerr << "n = " << n << "\n";
        return 4 * sqrt(lambda_b / n);

        // 2*sqrt(lambda_b/n) give semimajor axis. We want the length.

    }

    // return the new length that the cell would have
    // if site (x,y) were added.
    // used internally by CellularPotts
    inline double GetNewLengthIfXYWereAdded(int x, int y) const {
        double lengthifxywereadded = 0.;
        if (!par.periodic_boundaries)
            lengthifxywereadded = Length(sum_x + x, sum_y + y, sum_xx + x * x, sum_yy + y * y, sum_xy + x * y,
                                         area + 1);
        return lengthifxywereadded;
    }

    // return the new length that the cell would have
    // if site (x,y) were removed
    // used internally by CellularPotts
    inline double GetNewLengthIfXYWereRemoved(int x, int y) const {
        double lengthifxywereremoved = 0.;
        if (!par.periodic_boundaries)
            lengthifxywereremoved = Length(sum_x - x, sum_y - y, sum_xx - x * x, sum_yy - y * y, sum_xy - x * y,
                                           area - 1);
        return lengthifxywereremoved;
    }


private:
//! Increments the cell's actual area by 1 unit.
    inline int IncrementArea() {
        return ++area;
    }

    //! Decrements the cell's actual area by 1 unit.
    inline int DecrementArea() {
        return --area;
    }

    inline void DecrementAreaBy(int change) {
        area -= change;
    }

    /*! \brief Sets target area to actual area, to remove "pressure".

    This is useful when reading an initial condition from an image.
    */
    inline int SetAreaToTarget() {
        return area = target_area;
    }

    inline void SetTimeSinceBirth(int t) {
        time_since_birth = t;
    }

    inline int GetTimeSinceBirth() const {
        return time_since_birth;
    }

    // returns the maximum cell type index
    // (depends on Jtable)
    static int MaxTau() {
        return maxtau;
    }

protected:
    int group;
    bool alive;
    int sigma; // cell identity, 0 if medium
    int tau; // Cell type, when dynamicJ's are not used

    double meanx;
    double meany;
    //for Cell migration: target vector
    double tvecx;
    double tvecy;

    //stores old pos of cells for target vector calculations
    double prevx;
    double prevy;

    //save_movie direction of chemokine gradient (int plane)
    double chemvecx;
    double chemvecy;
    double chemmu; //this is the max strength of chemotaxis
    //migration parameters
    int persdur; //how long is this cell's persistent walk?
    int perstime; //counter for how long it has walked persistently
    double mu; //force of migration

    double length; // length of the cell;

    // so things have to be scaled a bit...

    //variables for dealing with genome and its output
    //mostly controlling when to divide right now
    Genome genome;
    int dividecounter; //how long have you had output saying to divide?
    int grad_conc; //how much gradient do I perceive?
    int gextiming; //funny variable deciding when to update gex



    // this is not used anymore - vJ is dynamically increased
    // Dynamically increased when cells are added to the system
    // unless a static Jtable is used (currently this is the default situation)
    static int **J;

    static int maxtau;

    // Amount: the number of Cell instantiations, INCLUDING copies
    // For internal use only.
    // Reading amount is NOT the way to get the number of cells!!
    static int amount;
    static int capacity;
    static int maxsigma; // the last cell identity number given out

    // TODO: Change to saving a vector of ancestors for better time resolution (see notes on notebook)
    // This is the oldest ancestor of this cell in the current season
    int ancestor;
    // indices of mother and daughter
    // (Note: no pointers, cells may be relocated)
    int mother;
    int daughter;
    int times_divided;
    int date_of_birth;
    int colour_of_birth;

    // These hold the decimal representation of the cell's jkey and jlock and is updated with the genome
    int jkey_dec = 0;
    int jlock_dec = 0;

    int area;
    int target_area;

    double v[2]{};
    int n_copies; // number of expansions of this cell
    // gradient of a chemical (to be extended to the total number chemicals)
    double grad[2]{};

    double *chem;
    // Raw moments of the cells
    // Are used to calculate minor and major axes
    // and center of mass
    // are locally adjusted, so axes are easily
    // and quickly calculated!

    // N.B: N is area!

    long int sum_x;
    long int sum_y;
    long int sum_xx;
    long int sum_yy;
    long int sum_xy;

    double food;
    int last_meal;

    int time_since_birth;

    const Dish *owner; // pointer to plane of cell

};

#endif
