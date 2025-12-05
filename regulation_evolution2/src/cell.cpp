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
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <map>

#ifndef __APPLE__

#include <malloc.h>

#endif

#include "cell.h"
#include "sticky.h"
#include "parameter.h"
#include "dish.h"

#define HASHCOLNUM 255

extern Parameter par;

int **Cell::J = nullptr;
int Cell::amount = 0;
int Cell::capacity = 0;
int Cell::maxsigma = 0;
int Cell::maxtau = 0;

//Cell::Cell(const Dish &who) : Cytoplasm(who);
// Note: g++ wants to have the body of this constructor in cell.hh
// body is defined in "ConstructorBody" below
class Dish;

using namespace std;

Cell::~Cell() {

    amount--;
    if (amount == 0) {
        // clear J if last cell has been destructed
        free(J[0]);
        free(J);
        capacity = 0;
        maxsigma = 0;
        J = nullptr;
    }
    if (par.n_chem) {
        delete[] chem;
    }
}


void Cell::CellBirth(Cell &mother_cell) {

    group = mother_cell.group;
    //cerr<< "daughter colour is "<<colour<<" mother colour is "<<mother_cell.colour<<endl;
    alive = mother_cell.alive;
    v[0] = mother_cell.v[0];
    v[1] = mother_cell.v[1];

    // Administrate ancestry
    ancestor = mother_cell.ancestor;
    mother_cell.daughter = this->sigma;
    mother = mother_cell.sigma;
    mother_cell.AddTimesDivided();
    times_divided = mother_cell.times_divided;

    owner = mother_cell.owner;

    date_of_birth = owner->Time();
    //cerr<<"sigma:"<<getSigma<<" I am born at: "<<date_of_birth<<endl<<endl; //this works fine

    colour_of_birth = mother_cell.group;
    group = mother_cell.group;

    alive = mother_cell.alive;

    tau = mother_cell.tau;
    target_area = mother_cell.target_area;

    genome = mother_cell.genome;
    gextiming = mother_cell.gextiming;
    dividecounter = 0;
    mother_cell.dividecounter = 0;
    //genome.MutateGenome();

    // Do not add moments here, they are going to be calculated from scratch
    // and are initialised to zero in ConstructorBody(), called right before this
    meanx = mother_cell.meanx;
    meany = mother_cell.meany;

    startTarVec();  //why this?
    mother_cell.startTarVec(); //randomises mother cell target vector
    // this fixes the biased movement bug!
    // also biologically it makes sense,
    // because cell polarity is all screwed after cell division
    //tvecx=mother_cell.tvecx;
    //tvecy=mother_cell.tvecy;

    prevx = mother_cell.prevx; //in startTarVec prevx is set to meanx, so we re-set this here
    prevy = mother_cell.prevy;
    persdur = mother_cell.persdur;
    perstime = int(persdur * RANDOM()); //assign random persistence duration upon birth.
    mu = mother_cell.mu;

    for (int ch = 0; ch < par.n_chem; ch++)
        chem[ch] = mother_cell.chem[ch];

    n_copies = 0;

    chemmu = mother_cell.chemmu;
    chemvecx = mother_cell.chemvecx;
    chemvecy = mother_cell.chemvecy;
    grad_conc = mother_cell.grad_conc;

    grad[0] = mother_cell.grad[0];
    grad[1] = mother_cell.grad[1];

    food = mother_cell.food;
    last_meal = mother_cell.last_meal;

    clearNeighbours(); //neighbours will be reassigned during the division function

    time_since_birth = 0;
    mother_cell.SetTimeSinceBirth(0);
}


void Cell::setNeighbour(int neighbour, int boundarylength, int contactduration) {

    if (boundarylength == 0)//remove this neighbour
        neighbours.erase(neighbour);
    else
        neighbours[neighbour] = make_pair(boundarylength,
                                          contactduration); //if the element is already present, the boundarylength will be modified, otherwise a new element will be created.

}

int Cell::returnBoundaryLength(int cell) {
    if (neighbours.count(cell))
        return neighbours[cell].first;

    return 0;
}


int Cell::returnDuration(int cell) {
    if (neighbours.count(cell))
        return neighbours[cell].second;

    return 0;
}

void Cell::clearNeighbours() {
    neighbours.clear();
}

int Cell::updateNeighbourBoundary(int cell, int boundarymodification) {
    //cerr<<"Hello updNeiBound begin"<<endl;

    if (!neighbours.count(cell) && boundarymodification < 0) {
        printf("Cell.updateNeighbourBoundary: error: negatively updating contact of cell %d with nonexisting neighbour %d\n",
               sigma, cell);
        exit(1);
        //return 1;
    } else if (!neighbours.count(cell)) {
        //cerr<<"Hello updNeiBound 0"<<endl;
        neighbours[cell] = make_pair(boundarymodification, 0);
        //cerr<<"Hello updNeiBound 1"<<endl;
    } else if (neighbours.count(cell)) {
        //cerr<<"Hello updNeiBound 2"<<endl;
        neighbours[cell].first += boundarymodification;
        //cerr<<"Hello updNeiBound 3"<<endl;
        //DO this after one MCS for duration reasons
//     if(neighbours[cell].first==0)//remove this neighbour
//     {
//       neighbours.erase(cell);
//     }
    }

    if (neighbours[cell].first < 0) {
        neighbours.erase(cell);
        printf("Cell.updateNeighbourBoundary: error: updating contact of cell %d with neighbour %d to negative value\n",
               sigma, cell);
        return 2;
    }

    return 0;
}

int Cell::SetNeighbourDurationFromMother(int cell, int motherduration) {
    if (!neighbours.count(cell)) {
        printf("Cell.SetNeighbourDurationFromMother: error: Nonexisting neighbour %d\n", cell);
        return 1;
    } else
        neighbours[cell].second = motherduration;

    return 0;
}


int Cell::updateNeighbourDuration(int cell, int durationmodification) {
    if (!neighbours.count(cell)) {
        printf("Cell.updateNeighbourDuration: error: Nonexisting neighbour %d\n", cell);
        return 1;
    } else if (neighbours.count(cell)) {
        neighbours[cell].second += durationmodification;
    }


    return 0;
}

// TODO change this to an iterator that uses a spiral to get cell positions
BoundingBox Cell::getBoundingBox() {
    double maxx = owner->CPM->SizeX() - 1;
    double maxy = owner->CPM->SizeY() - 1;
    // Veryyy arbitrary scaling to make sure all positions are in
    double scale = 3;
    return BoundingBox{
            (int) max(meanx - scale * Length(), 0.),
            (int) max(meany - scale * Length(), 0.),
            (int) min(meanx + scale * Length(), maxx),
            (int) min(meany + scale * Length(), maxy)
    };
}


void Cell::updateJDecs() {
    auto p = genome.calculateJdecs();
    jkey_dec = p.first;
    jlock_dec = p.second;
}
