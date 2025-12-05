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

#ifndef __APPLE__

#include <malloc.h>

#endif

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "dish.h"
#include "random.h"
#include "output.h"
#include "misc.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else

#include "x11graph.h"

#endif

//NOTE: The bookkeeping for cell contacts is very extensive:
// When cells are initially placed (call dish->InitContactLength afterwards)
// When cells divide (in cpm->dividecells)
// When cells are killed and removed (cpm->removecells)
// During CPM updates (cpm->convertspin)
// I added pieces of code to take care of this in the various applicable functions
// We may want to add a parameter to make these parts optional in case we don't need it --it's a bit more costly

using namespace std;

INIT {
    try {
        if (strlen(par.celldatafile) and strlen(par.latticefile)) {
            cout << "Reading backfile" << endl;
            cout << "backup file is " << par.latticefile << endl;
            setStartTime(readCellData());
            cout << par.fooddatafile << endl;
            if (!strlen(par.fooddatafile)) {
                for (int i = 0; i < par.foodpatches; i++) {
                    food_manager.addRandomFoodPatch();
                }
            } else if (readFoodData() != getStartTime())
                cerr << "Food data and cell data date from different times!" << endl;
            readLattice();
            CPM->InitializeEdgeList(false);
            InitContactLength();
        } else {
            //THIS IS TO USE FOR NORMAL INITIALISATION
            if (par.scatter_start) {
                CPM->PlaceCellsRandomly(par.n_init_cells, par.size_init_cells);
            } else {
                CPM->PlaceCellsOrderly(par.n_init_cells, par.size_init_cells);
            }
            CPM->InitializeEdgeList(false);
            cout << "done initialising edge list" << endl;

            CPM->ConstructInitCells(*this); //within an object, 'this' is the object itself

            cout << "done setting types" << endl;
            //Initialise the contactlength bookkeeping now that the cells are placed
            // at this stage, cells are only surrounded by medium
            InitContactLength();  // see dish.cpp - you don't need dish->InitContactLength because this part IS in dish
            cout << "done setting contact length" << endl;

            cout << "Going to initialise genome" << endl;
            for (auto &c: cell) {
                if (c.Sigma()) {
                    c.setGTiming((int) (RANDOM() * par.grn_update_period));
                    //c.setGTiming(0);
                    c.dividecounter = 0;
                    c.SetTargetArea(
                        par.target_area); //sets target area because in dividecells the new target area = area
                    //initialise a cell's timing for gex Updating
                    //c.setGTiming((int)(RANDOM()*par.grn_update_period));
                    //creates a cell's genome, either randomly or from file
                    if (strlen(par.genomefile)) {
                        c.ReadGenomeFromFile(par.genomefile);
                    } else {
                        int key_lock_len = int(stringToVector<int>(par.key_lock_weights, ' ').size());
                        c.CreateRandomGenome(
                                2,
                                par.nr_regnodes,
                                1 + key_lock_len * 2,
                                stringToVector<double>(par.input_scales, ' ')
                        );
                    }
                    c.ClearGenomeState();
                }
            }

            cout << "Done initialising genome" << endl;
            for (int i = 0; i < par.foodpatches; ++i)
                food_manager.addRandomFoodPatch();
            cout << "done with food" << endl;
            //run CPM for some time without persistent motion
            for (int init_time = 0; init_time < 10; init_time++) {
                CPM->AmoebaeMove2(PDEfield);  //this changes neighs
            }
            cout << "done with update" << endl;
            InitCellMigration();
            UpdateCellParameters(0);//update cell status //UpdateCellParameters2();
            setStartTime(0);
        }
    } catch (const char *error) {
        cerr << "Caught exception\n";
        std::cerr << error << "\n";
        exit(1);
    }
}


TIMESTEP {
    try {
        static Dish *dish = new Dish(); //here ca planes and cells are constructed
        static int i = dish->getStartTime(); //starttime is set in Dish. Not the prettiest solution, but let's hope it works.

        dish->cellsEat(i);

        dish->UpdateCellParameters(i); // for continuous GRN updating and reproduction

        dish->CellMigration();//updates persistence time and targetvectors

        dish->CPM->AmoebaeMove2(dish->PDEfield);  //this changes neighs

        dish->UpdateNeighDuration();

        if (i % 25 == 0) {
            // Try to add fpatch
            dish->getFood().replenishFood(i);

            //check if one of the groups is extinct. if yes, end simulation
            if (par.groupextinction and dish->groupExtinction()) {
                std::cout << "Group extinct after " << i << " time steps. ending simulation...\n";
                return false;
            }
        }

        if (i % par.info_period == 0) {
            cout << "Time = " << i << '\n';
            cout << "There are " << dish->CountCells() << " cells" << endl;  // Time to flush
        }

        // TO FILE FOR MOVIE
        if (par.save_movie && !(i % par.save_movie_period)) {
            dish->makePlots(i, this);
        }
        if (!(i % par.save_data_period)) {
            dish->saveFoodData(i);
            if (not dish->getCellGraves().empty())
                dish->saveCellGraveData(i);
            int popsize = dish->saveCellData(i);
            if (not popsize) {
                cout << "Global extinction after " << i << " time steps, simulation terminates now" << endl;
                return false;
            }
        }
        // TO FILE FOR BACKUP
        if (!(i % par.save_lattice_period)) {
            dish->saveLattice(i);
        }

        i++;
        return i <= par.mcs;
    } catch (const char *error) {
        cerr << "Caught exception\n";
        std::cerr << error << "\n";
        exit(1);
    }
}

int PDE::MapColour(double val) {

    return (((int) ((val / ((val) + 1.)) * 100)) % 100) + 155;
}

//////////////////////////////
// ------------------------ //
// ---       MAIN       --- //
// ------------------------ //
//////////////////////////////
int main(int argc, char *argv[]) {
    try {
        par.Read(argv[1]); // Read parameters from file

        //command line arguments overwrite whatever is in the parameter file
        if (argc > 2) {
            int exit_valarg = par.ReadArguments(argc, argv);
            if (0 != exit_valarg) {
                Parameter::PrintWelcomeStatement(); //see parameter.h
                return EXIT_FAILURE;
            }
        }

        cerr << endl << "Warning, this version is ***NOT*** suitable for pde field!!!" << endl;
        //Depends on this: AddSiteToMoments (and Remove), FindCellDirections2, etc...
        cerr << endl << "WARNING, use wrapped boundaries if cells are A LOT smaller than sizex and sizey" << endl
             << endl;
        cerr << endl
             << "WARNING: DO NOT EVOLVE CHEMMU, or if you do, change the replication function (where it is always reset to init_chemmu)"
             << endl << endl;

        if (par.existing_dirs) {
            cerr << "WARNING: Outputting data to existing directory structure\n";
            if (par.replace_dirs)
                cerr << "WARNING: Deleting previously existing data\n";
        }
        //check if directory for movies exists, create it if not, exit otherwise
        makeDirIfNeeded(par.moviedir, par.existing_dirs, par.replace_dirs);  //see output.cpp
        makeDirIfNeeded(par.latticedir, par.existing_dirs, par.replace_dirs);  //see output.cpp
        makeDirIfNeeded(par.celldatadir, par.existing_dirs, par.replace_dirs);  //see output.cpp
        makeDirIfNeeded(par.cellgravesdatadir, par.existing_dirs, par.replace_dirs);  //see output.cpp
        makeDirIfNeeded(par.fooddatadir, par.existing_dirs, par.replace_dirs);  //see output.cpp

        std::cout << "Seed for the random generator is: " << Seed(par.rseed) << "\n";

        cout << "Using X11 graphics (batch mode). sizex and y are " << par.sizex << " " << par.sizey << endl;
        X11Graphics g(par.sizex, par.sizey);

        // Runs at most par.mcs steps (can stop early)
        bool running = true;
        while (running) {
            running = g.TimeStep();
        }
        cout << "End of simulation, goodbye!";
        return 0;
    } catch (const char *error) {
        std::cerr << error << "\n";
        return 1;
    }
}
