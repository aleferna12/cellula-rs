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
#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

class Parameter {
public:
    Parameter();

    ~Parameter();

    void CleanUp(void) const;

    int ReadArguments(int argc, char *argv[]);

    static void PrintWelcomeStatement(void);

    void Read(const char *filename);

    void Write(ostream &os) const;

    // Boltzmann temperature.
    double T = 16.;
    // Target area constraint in the hamiltonian.
    int target_area = 50;
    // Size lambda.
    double lambda = 4;
    // Whether to penalize loose pixels.
    int conn_diss = 2000;
    // Energy of cell interactions with border.
    int border_energy = 100;
    // Degree of neighbourhood interactions (2 = Moore neighs).
    int neighbours = 2;
    // Minimum area allowed for a living cell.
    int min_area_for_life = 5;
    // Weights of cell-cell adhesion.
    string key_lock_weights = "1 2 3 4 5 6";
    // Adhesion constant for cell-medium interactions.
    double Jmed = 14;
    // Adhesion constant for cell-cell interactions.
    double Jalpha = 7;
    int n_chem = 0;
    double *diff_coeff;
    double *decay_rate;
    double *secr_rate;
    double dt = 2.0;
    double dx = 2.0e-6;
    // Number of cells initialized in the beginning of a sim.
    int n_init_cells = 100;
    // Size of the cells initialized in the beginning of a sim.
    int size_init_cells = 25;
    // X component of the lattice size (only square lattices are allowed for now).
    int sizex = 200;
    // Y component of the lattice size (only square lattices are allowed for now).
    int sizey = 200;
    // Number of MCSs to run the simulation for before it's interrupted.
    int mcs = 10000;
    // Random generator seed (positive values are seeds, -1 means "pick random seed").
    unsigned int rseed = -1;
    // Whether to interrupt the sim in case one of the groups gets extinct in competition sims.
    bool groupextinction = false;
    // Whether cell groups are initialized randomly at the start of the sim. By default, all cells' group is '0'.
    bool randgroups = false;
    // Whether adhesion is determined dynamically by the cells' GRNs.
    // If false, the Jcell-cell term of the adhesion gamma equation simplifies to "Jalpha".
    bool dynamic_adh = true;
    // File containing a genome that will be initialized in the starting cells of the sim.
    char *genomefile = strdup("");
    // How many regulatory nodes in genomes (should be around number of output nodes + 1).
    int nr_regnodes = 1;
    // How many steps it takes before a cell can start dividing.
    int divtime = 20;
    // Duration of division once a cell starts division program.
    int divdur = 1000;
    // Chance of the genome being mutated.
    double mu = 0.05;
    // Standard deviation for the mutation rate of the genome.
    double mustd = 0.1;
    // Alpha term of the Gompertz's hazard function used to kill cells.
    double gompertz_alpha = 0.0000075;
    // Beta term of the Gompertz's hazard function used to kill cells.
    // If set to zero cells are killed independently of age with a random probability determined by "gompertz_alpha".
    double gompertz_beta = 0.0004;
    // How much food do cells start with.
    int initial_food_amount = 0;
    // How often in MCSs do we take away one food particle from the cells.
    int metabperiod = 20;
    // Noise associated with the edges of the chemotactic gradient
    double gradnoise = 0.1;
    // Controls the depth of the gradient (smaller "gradscale", shallower gradient).
    // gradscale = 100 / gradient_step_size
    double gradscale = 1.0;
    // Number of food patches initially spawned on the field.
    int foodpatches = 1;
    // How often in MCSs we spawn a new food patch (as long as we have less than "maxfoodpatches").
    // Seasonality might affect this rate.
    int foodpatchperiod = 1000;
    // How long a season cycle lasts for (period of the sine wave).
    int seasonduration = 10000;
    // What is the amplitude of the seasonal variation (amplitude of the sine wave).
    double seasonamplitude = 0;
    // Diameter of each food patch (total amount of spots per patch will be approx. = pi * (length/2)^2).
    int foodpatcharea = 1;
    // How much food each spot in the food patch contains.
    int foodperspot = 1;
    // How much food does each cell start with.
    int foodstart = 1000;
    // How often in MCSs are cells allowed to eat 1 food.
    int eatperiod = 1000;
    // Maximum number of food patches allowed to coexist in the system.
    int maxfoodpatches = 1;
    // Whether to not allow cell division.
    bool nodivisions = false;
    // Whether to not allow cell death.
    bool nodeaths = false;
    // Chance that cells copy medium from nowhere instead of from each other (prevents extinction of medium cell).
    double chancemediumcopied = 0.0001;
    // Relative path to the colo table file.
    // Use the scripts/colortable module to generate these from colorir palettes.
    string colortablefile = "default.ctb";
    // Which plots to make, specified as a space-delimited string.
    // Possible values for the list are: tau food group.
    string plots = "tau";
    // Inputs are multiplied by these values before being passed to the network
    // You should pick values that make the order of magnitude of the inputs match
    string input_scales = "0.1 1";
    // Whether to draw the lines showing the direction of migration for each cell
    bool miglines = true;
    // Whether to draw the chemotactic gradient.
    bool chemgrad = true;
    // Whether to draw concentric circles around the food.
    bool chemcircles = true;
    // Thickness of the gradient circle lines, in pixels.
    int circle_thickness = 1;
    // Distance between the gradient circle lines, in pixels.
    int circle_dist = 40;
    // Number of black line segments drawn for the most proximal gradient line (second will have double that amount
    // and so on).
    int circle_segments = 15;
    // How much of the circumference of the gradient circle lines is black (interval [0, 1]).
    double circle_coverage = 0.5;
    double startmu = 0.0;
    double init_chemmu = 0.;
    int persduration = 0;
    // How often in MCSs do we update the GRN of cells (and also any properties that depend on that such as tau).
    int grn_update_period = 1;
    // Whether the initial pop. of cells should be scattered in the beginning of the sim. or rather start together.
    bool scatter_start = true;
    // Relative path to the directory where lattice backup data will be saved.
    char *latticedir = strdup("lattice");
    // Relative path to the directory where cell data will be saved.
    char *celldatadir = strdup("celldata");
    // Relative path to the directory where food data will be saved.
    char *fooddatadir = strdup("fooddata");
    // Relative path to the directory where cell death data will be saved.
    char *cellgravesdatadir = strdup("cellgravedata");
    // Where to store sim pictures.
    char *moviedir = strdup("data_film");
    // Relative path to the lattice file used to reinitialize the sim. from a backup.
    char *latticefile = strdup("");
    // Relative path to the cell data file used to reinitialize the sim. from a backup.
    char *celldatafile = strdup("");
    // Relative path to the food data file used to reinitialize the sim. from a backup.
    char *fooddatafile = strdup("");
    // Whether it is ok to output data to directories that already exist
    // Can lead to existing data being overwritten
    bool existing_dirs = false;
    // Whether to delete existing data directories before starting the simulation
    // Only works if 'existing_dirs' is true
    bool replace_dirs = false;
    // How often in MCSs to save cell, cell death and food data.
    int save_data_period = 100;
    // How often in MCSs to save lattice backup files.
    int save_lattice_period = 1000;
    // How often to save pictures of the lattice.
    int save_movie_period = 100;
    // Whether to save pictures of the lattice.
    bool save_movie = false;
    // How often in MCSs to output information about the sim. in the terminal.
    int info_period = 1000;
    // Whether to plot simulation in real time (currently not supported).
    bool graphics = true;
    // Whether to mutate the GRN.
    bool evolreg = true;
    // Initial configuration of cells, see ca.cpp placecellsorderly().
    int init_cell_config = 0;
    int cell_placement = 0;
    // Whether to use periodic boundaries to wrap the sim.
    // No longer supported (should always be false) but kept in case we want to come back to it at some point.
    bool periodic_boundaries = false;
};

ostream &operator<<(ostream &os, Parameter &p);

const char *sbool(const bool &p);

extern Parameter par;

#endif
