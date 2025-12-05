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


#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "output.h"
#include "parse.h"

Parameter::Parameter() {
    diff_coeff = new double[1];
    diff_coeff[0] = 1e-13;
    decay_rate = new double[1];
    decay_rate[0] = 1.8e-4;
    secr_rate = new double[1];
}

Parameter::~Parameter() {

    // destruct parameter object

    // free string parameter

    CleanUp();

}

void Parameter::CleanUp() const {
    if (diff_coeff)
        free(diff_coeff);
    if (decay_rate)
        free(decay_rate);
    if (secr_rate)
        free(secr_rate);

    if (moviedir)
        free(moviedir);
    if (genomefile)
        free(genomefile);
    if (latticedir)
        free(latticedir);
    if (celldatadir)
        free(celldatadir);
    if (cellgravesdatadir)
        free(cellgravesdatadir);
    if (fooddatadir)
        free(fooddatadir);
    if (latticefile)
        free(latticefile);
}

void Parameter::PrintWelcomeStatement() {
    cout << "CellEvol: v0.something (very much a prototype)" << endl;
    cout << "Usage is: " << endl;
    cout << "./cell_evolution path/to/data [optional arguments]" << endl;
    cout << "Arguments: " << endl;
    cout
            << " -name path/to/name_for_all_output # gives a name to all output, alternative to -moviedir -latticedir -celldatadir -cellgravedatadir -fooddatadir"
            << endl;
    cout << " -celldatafile path/to/celldatafile # output file" << endl;
    cout << " -fooddatafile path/to/fooddatafile # output file" << endl;
    cout << " -moviedir path/to/moviedir # output movie dir" << endl;
    cout << " -latticedir path/to/latticedir # output backup lattice dir" << endl;
    cout << " -celldatadir path/to/celldatadir # output cell data dir" << endl;
    cout << " -cellgravedatadir path/to/cellgravedatadir # output cell data dir" << endl;
    cout << " -fooddatadir path/to/fooddatadir # output food data dir" << endl;
    cout << " -save_movie # save_movie pictures" << endl;
    cout << " -seed INT_NUMBER # for random number generator" << endl;
    cout << " -mcs INT_NUMBER" << endl;
    // cout<<" -halfdiv_area_predator INT_NUMBER"<<endl;
    cout << " -persmu FLOAT_NUMBER [ > 0 ], strength of persistent random walk" << endl;
    cout << " -persduration INT_NUMBER" << endl;
    cout << " -mu FLOAT_NUMBER [0,1) # mutation rate for genome regulation" << endl;
    cout << " -mustd FLOAT_NUMBER [0,1) # mutation size for genome regulation" << endl;
    cout << " -casize INT_NUMBER INT_NUMBER # dimensions of the CA" << endl;
    cout << " -noevolreg # No evolution of regulation parameters" << endl;
    cout << " -scatter # spread cells after a season" << endl;
    cout << " -nodivisions # do not execute divisions -> number of cells remains the same" << endl;
    cout << " -nodeaths # do not kill cells -> number of cells remains the same" << endl;
    cout << " -latticefile path/to/latticefile # to start simulation from backup" << endl;
    cout << " -colortablefile path/to/colortablefile # to specify color table" << endl;
    cout << " -info_period [INT_NUMBER] # how often output sim information" << endl;
    cout << " -metabperiod [INT_NUMBER] how often we deduce 1 food of each cell in MCS" << endl;
    cout << " -gradscale [FLOAT_NUMBER] slope of the gradient (in percent units)" << endl;
    cout << " -foodpatches [INT_NUMBER] initial number of food patches resources placed in the field" << endl;
    cout << " -foodpatchperiod [INT_NUMBER] new food patch timer (a new patch will be created every X MCS)" << endl;
    cout << " -seasonduration [INT_NUMBER] length of the season cycle" << endl;
    cout << " -seasonamplitude [FLOAT_NUMBER] amplitude of the seasonal variation" << endl;
    cout << " -foodpatcharea [INT_NUMBER] are of each food patch" << endl;
    cout << " -foodperspot [INT_NUMBER] how much food each spot contains" << endl;
    cout << " -maxfoodpatches [INT_NUMBER] maximum number of food patches allowed to coexist in the system" << endl;
    cout << " -foodstart [INT_NUMBER] initial food for cells" << endl;
    cout << " -eatperiod [INT_NUMBER] how often a cell can eat" << endl;
    cout << " -gradnoise [FLOAT_NUMBER] chances that any grid point has gradient, rather than being empty" << endl;
    cout << " -chemmu [FLOAT_NUMBER] scaling factor for chemotaxis in the Hamiltonian" << endl;
    cout << " -genomefile [string] starting genome with which to seed the field" << endl;
    cout << " -target_area [INT_NUMBER] that (initial) target area of cells" << endl;
    cout << " -init_cell_config [0-3] initial configuration of cells when placed in center, see ca.cpp" << endl;
    cout << " -cell_placement [1-4] field position of cells, (0=center) see ca.cpp" << endl;
    cout << " -gompertz_alpha [FLOAT_NUMBER] alpha parameter for the hazard function a*e^(b*x)" << endl;
    cout << " -gompertz_beta [FLOAT_NUMBER] beta parameter for the hazard function a*e^(b*x)" << endl;
    cout << endl << "Will not execute if celldatafile and moviedir already exist" << endl;
    cout << "Also, parameter file and Jtable should be in the same directory (unless you used option -keylockfilename)"
         << endl;
    cout << "Have fun!" << endl;
}

int Parameter::ReadArguments(int argc, char *argv[]) {
    cout << endl << "Reading arguments from command line" << endl;
    //starts from 2 because 0 is filename, 1 is parameter file path
    for (int i = 2; i < argc; i++) {
        if (0 == strcmp(argv[i], "-celldatafile")) {
            i++;
            if (i == argc) {
                cout << "Something odd in celldatafile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(celldatafile, argv[i]); //this can be buggy because it copies over a dynamically allocated char* (celldatafile) that can be a lot shorter
            free(celldatafile);
            celldatafile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            celldatafile = strdup(argv[i]);
            cout << "New value for celldatafile: " << celldatafile << endl;
//       exit(1);
        } else if (0 == strcmp(argv[i], "-fooddatafile")) {
            i++;
            if (i == argc) {
                cout << "Something odd in fooddatafile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(fooddatafile, argv[i]); //this can be buggy because it copies over a dynamically allocated char* (celldatafile) that can be a lot shorter
            free(fooddatafile);
            fooddatafile = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            fooddatafile = strdup(argv[i]);
            cout << "New value for fooddatafile: " << fooddatafile << endl;
        } else if (0 == strcmp(argv[i], "-moviedir")) {
            i++;
            if (i == argc) {
                cout << "Something odd in moviedir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(moviedir);
            moviedir = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            moviedir = strdup(argv[i]);

            cout << "New value for moviedir: " << moviedir << endl;

        } else if (0 == strcmp(argv[i], "-latticedir")) {
            i++;
            if (i == argc) {
                cout << "Something odd in latticedir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(latticedir);
            latticedir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            latticedir = strdup(argv[i]);

            cout << "New value for latticedir: " << latticedir << endl;

        } else if (0 == strcmp(argv[i], "-celldatadir")) {
            i++;
            if (i == argc) {
                cout << "Something odd in celldatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(celldatadir);
            celldatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            celldatadir = strdup(argv[i]);

            cout << "New value for celldatadir: " << celldatadir << endl;

        } else if (0 == strcmp(argv[i], "-cellgravesdatadir")) {
            i++;
            if (i == argc) {
                cout << "Something odd in cellgravesdatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(cellgravesdatadir);
            cellgravesdatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            cellgravesdatadir = strdup(argv[i]);

            cout << "New value for cellgravesdatadir: " << cellgravesdatadir << endl;

        } else if (0 == strcmp(argv[i], "-fooddatadir")) {
            i++;
            if (i == argc) {
                cout << "Something odd in fooddatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(fooddatadir);
            fooddatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            fooddatadir = strdup(argv[i]);

            cout << "New value for fooddatadir: " << fooddatadir << endl;
        } else if (0 == strcmp(argv[i], "-genomefile")) {
            i++;
            if (i == argc) {
                cout << "Something odd in genomefile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(genomefile);
            genomefile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            genomefile = strdup(argv[i]);

            cout << "New value for genomefile: " << genomefile << endl;

        } else if (0 == strcmp(argv[i], "-seed")) {
            i++;
            if (i == argc) {
                cout << "Something odd in seed?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            rseed = atoi(argv[i]);
            cout << "New value for seed: " << rseed << endl;
        } else if (0 == strcmp(argv[i], "-mcs")) {
            i++;
            if (i == argc) {
                cout << "Something odd in mcs?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mcs = atoi(argv[i]);
            cout << "New value for mcs: " << mcs << endl;
        } else if (0 == strcmp(argv[i], "-circle_segments")) {
            i++;
            if (i == argc) {
                cout << "Something odd in circle_segments?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_segments = atoi(argv[i]);
            cout << "New value for circle_segments (mcs in the code): " << circle_segments << endl;
        } else if (0 == strcmp(argv[i], "-circle_thickness")) {
            i++;
            if (i == argc) {
                cout << "Something odd in circle_thickness?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_thickness = atoi(argv[i]);
            cout << "New value for circle_thickness (mcs in the code): " << circle_thickness << endl;
        } else if (0 == strcmp(argv[i], "-circle_dist")) {
            i++;
            if (i == argc) {
                cout << "Something odd in circle_dist?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_dist = atoi(argv[i]);
            cout << "New value for circle_dist (mcs in the code): " << circle_dist << endl;
        } else if (0 == strcmp(argv[i], "-circle_thickness")) {
            i++;
            if (i == argc) {
                cout << "Something odd in circle_thickness?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_thickness = atof(argv[i]);
            cout << "New value for circle_thickness: " << circle_thickness << endl;
        } else if (0 == strcmp(argv[i], "-mu")) {
            i++;
            if (i == argc) {
                cout << "Something odd in mu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mu = atof(argv[i]);
            cout << "New value for genome mutation rate: " << mu << endl;
        } else if (0 == strcmp(argv[i], "-mustd")) {
            i++;
            if (i == argc) {
                cout << "Something odd in mustd?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mustd = atof(argv[i]);
            cout << "New value for mustd: " << mustd << endl;
        } else if (0 == strcmp(argv[i], "-persduration")) {
            i++;
            if (i == argc) {
                cout << "Something odd in persduration?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            persduration = atoi(argv[i]);
            cout << "New value for persistence of movement: " << persduration << endl;
        } else if (0 == strcmp(argv[i], "-latticefile")) {
            i++;
            if (i == argc) {
                cout << "Something odd in latticefile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            free(latticefile);
            latticefile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            latticefile = strdup(argv[i]);

            cout << "New value for latticefile: " << latticefile << endl;
        } else if (0 == strcmp(argv[i], "-colortablefile")) {
            i++;
            if (i == argc) {
                cout << "Something odd in colortablefile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            colortablefile = argv[i];

            cout << "New value for colortablefile: " << colortablefile << endl;
        } else if (0 == strcmp(argv[i], "-casize")) {
            i++;
            if (i == argc) {
                cout << "Something odd in casize?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            sizex = atoi(argv[i]);
            i++;
            if (i == argc) {
                cout << "Something odd in casize?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            sizey = atoi(argv[i]);
            cout << "New value for CA size x and y: " << sizex << " " << sizey << endl;
        } else if (0 == strcmp(argv[i], "-info_period")) {
            i++;
            if (i == argc) {
                cout << "Something odd in info_period?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            info_period = atoi(argv[i]);
            cout << "New value for info_period: " << info_period << endl;
        } else if (0 == strcmp(argv[i], "-noscatter_start")) {
            scatter_start = false;
            cout << "Cells will not be scattered at the start of the first season" << endl;
        } else if (0 == strcmp(argv[i], "-nochemgrad")) {
            chemgrad = false;
            cout << "Not plotting chemotactic gradient" << endl;
        } else if (0 == strcmp(argv[i], "-groupextinction")) {
            groupextinction = true;
            cout << "Stopping sim on group extinction" << endl;
        } else if (0 == strcmp(argv[i], "-randgroups")) {
            randgroups = true;
            cout << "Random group initialization" << endl;
        }  else if (0 == strcmp(argv[i], "-noevolvable_adh")) {
            dynamic_adh = false;
            cout << "Adhesion is not evolvable" << endl;
        } else if (0 == strcmp(argv[i], "-nochemcircles")) {
            chemcircles = false;
            cout << "Not plotting circles around food patches" << endl;
        } else if (0 == strcmp(argv[i], "-nomiglines")) {
            miglines = false;
            cout << "Not plotting migration direction lines" << endl;
        } else if (0 == strcmp(argv[i], "-existing_dirs")) {
            existing_dirs = true;
            cout << "Not plotting circles around food patches" << endl;
        } else if (0 == strcmp(argv[i], "-replace_dirs")) {
            replace_dirs = true;
            cout << "Not plotting circles around food patches" << endl;
        } else if (0 == strcmp(argv[i], "-nodivisions")) {
            nodivisions = true;
            cout << "Cells will not ever die" << endl;
        } else if (0 == strcmp(argv[i], "-nodeaths")) {
            nodeaths = true;
            cout << "Cells will not ever die" << endl;
        } else if (0 == strcmp(argv[i], "-save_movie")) {
            save_movie = true;
            cout << "pictures will be stored" << endl;
        } else if (0 == strcmp(argv[i], "-noevolreg")) {
            evolreg = false;
            cout << "No evolution of regulation parameters" << endl;
        } else if (0 == strcmp(argv[i], "-metabperiod")) {
            i++;
            if (i == argc) {
                cout << "Something odd in metabperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            metabperiod = atoi(argv[i]);
            cout << "New value for metabperiod: " << metabperiod << endl;
        } else if (0 == strcmp(argv[i], "-save_movie_period")) {
            i++;
            if (i == argc) {
                cout << "Something odd in save_movie_period?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            save_movie_period = atoi(argv[i]);
            cout << "New value for save_movie_period: " << save_movie_period << endl;
        } else if (0 == strcmp(argv[i], "-save_lattice_period")) {
            i++;
            if (i == argc) {
                cout << "Something odd in save_lattice_period?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            save_lattice_period = atoi(argv[i]);
            cout << "New value for save_lattice_period: " << save_lattice_period << endl;
        } else if (0 == strcmp(argv[i], "-save_data_period")) {
            i++;
            if (i == argc) {
                cout << "Something odd in save_data_period?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            save_data_period = atoi(argv[i]);
            cout << "New value for save_data_period: " << save_data_period << endl;
        } else if (0 == strcmp(argv[i], "-gradscale")) {
            i++;
            if (i == argc) {
                cout << "Something odd in gradscale?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            gradscale = atof(argv[i]);
            cout << "New value for gradscale: " << gradscale << endl;
        } else if (0 == strcmp(argv[i], "-foodpatches")) {
            i++;
            if (i == argc) {
                cout << "Something odd in foodpatches?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatches = atoi(argv[i]);
            cout << "New value for foodpatches: " << foodpatches << endl;
        } else if (0 == strcmp(argv[i], "-foodpatchperiod")) {
            i++;
            if (i == argc) {
                cout << "Something odd in foodpatchperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatchperiod = atoi(argv[i]);
            cout << "New value for foodpatchperiod: " << foodpatchperiod << endl;
        } else if (0 == strcmp(argv[i], "-seasonduration")) {
            i++;
            if (i == argc) {
                cout << "Something odd in seasonduration?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            seasonduration = atoi(argv[i]);
            cout << "New value for seasonduration: " << seasonduration << endl;
        } else if (0 == strcmp(argv[i], "-seasonamplitude")) {
            i++;
            if (i == argc) {
                cout << "Something odd in seasonamplitude?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            seasonamplitude = atof(argv[i]);
            cout << "New value for seasonamplitude: " << seasonamplitude << endl;
        } else if (0 == strcmp(argv[i], "-foodpatcharea")) {
            i++;
            if (i == argc) {
                cout << "Something odd in foodpatcharea?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatcharea = atoi(argv[i]);
            cout << "New value for foodpatcharea: " << foodpatcharea << endl;
        } else if (0 == strcmp(argv[i], "-foodperspot")) {
            i++;
            if (i == argc) {
                cout << "Something odd in foodperspot?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodperspot = atoi(argv[i]);
            cout << "New value for foodperspot: " << foodperspot << endl;
        } else if (0 == strcmp(argv[i], "-maxfoodpatches")) {
            i++;
            if (i == argc) {
                cout << "Something odd in maxfoodpatches?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            maxfoodpatches = atoi(argv[i]);
            cout << "New value for maxfoodpatches: " << maxfoodpatches << endl;
        } else if (0 == strcmp(argv[i], "-foodstart")) {
            i++;
            if (i == argc) {
                cout << "Something odd in foodstart?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodstart = atoi(argv[i]);
            cout << "New value for foodstart: " << foodstart << endl;
        } else if (0 == strcmp(argv[i], "-eatperiod")) {
            i++;
            if (i == argc) {
                cout << "Something odd in eatperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            eatperiod = atoi(argv[i]);
            cout << "New value for eatperiod: " << eatperiod << endl;
        } else if (0 == strcmp(argv[i], "-chemmu")) {
            i++;
            if (i == argc) {
                cout << "Something odd in chemmu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            init_chemmu = atof(argv[i]);
            cout << "New value for chemmu: " << init_chemmu << endl;
        } else if (0 == strcmp(argv[i], "-target_area")) {
            i++;
            if (i == argc) {
                cout << "Something odd in target_area?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            target_area = atoi(argv[i]);
            cout << "New value for target_area: " << target_area << endl;
        } else if (0 == strcmp(argv[i], "-persmu")) {
            i++;
            if (i == argc) {
                cout << "Something odd in persmu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            startmu = atof(argv[i]);
            cout << "New value for persmu: " << startmu << endl;
        } else if (0 == strcmp(argv[i], "-init_cell_config")) {
            i++;
            if (i == argc) {
                cout << "Something odd in init_cell_config?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            init_cell_config = atoi(argv[i]);
            cout << "New value for init_cell_config: " << init_cell_config << endl;
        } else if (0 == strcmp(argv[i], "-cell_placement")) {
            i++;
            if (i == argc) {
                cout << "Something odd in cell_placement?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            cell_placement = atoi(argv[i]);
            cout << "New value for cell_placement: " << cell_placement << endl;
        } else if (0 == strcmp(argv[i], "-gradnoise")) {
            i++;
            if (i == argc) {
                cout << "Something odd in gradnoise?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            gradnoise = atof(argv[i]);
            cout << "New value for gradnoise: " << gradnoise << endl;
        } else if (0 == strcmp(argv[i], "-name")) {
            i++;
            if (i == argc) {
                cout << "Something odd in name?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            // I'm just going to work in c++ strings - a lot easier
            free(moviedir);
            free(latticedir);
            free(celldatadir);
            free(cellgravesdatadir);
            free(fooddatadir);

            string maybepath_and_name(argv[i]);
            size_t botDirPos = maybepath_and_name.find_last_of("/");
            string dir("");
            string name;
            if (botDirPos != std::string::npos) {
                // then there is a character '/' in name, which means that
                // we are going to save data in some path, hence
                // we have to split where this is happening
                dir = maybepath_and_name.substr(0, botDirPos + 1);
                name = maybepath_and_name.substr(botDirPos + 1, maybepath_and_name.length());
            } else {
                name = maybepath_and_name;
            }

            string name_outfile = dir; //will this contain the last '/''
            name_outfile.append("data_");
            name_outfile.append(name);
            name_outfile.append(".txt");

            string name_moviedir = dir;
            name_moviedir.append("movie_");
            name_moviedir.append(name);

            string name_latticedir = dir;
            name_latticedir.append("lattice_");
            name_latticedir.append(name);

            string name_celldatadir = dir;
            name_celldatadir.append("celldata_");
            name_celldatadir.append(name);

            string name_cellgravedatadir = dir;
            name_cellgravedatadir.append("cellgravedata_");
            name_cellgravedatadir.append(name);

            string name_fooddatadir = dir;
            name_fooddatadir.append("fooddata_");
            name_fooddatadir.append(name);

            std::cout << "New value for output directory: " << name_outfile << '\n';
            std::cout << "New value for moviedir: " << name_moviedir << '\n';
            std::cout << "New value for latticedir: " << name_latticedir << '\n';
            std::cout << "New value for celldatadir: " << name_celldatadir << '\n';
            std::cout << "New value for cellgravedatadir: " << name_cellgravedatadir << '\n';
            std::cout << "New value for fooddatadir: " << name_fooddatadir << '\n';

            moviedir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            latticedir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            celldatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            cellgravesdatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            fooddatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            strdup(name_outfile.c_str());
            moviedir = strdup(name_moviedir.c_str());
            latticedir = strdup(name_latticedir.c_str());
            celldatadir = strdup(name_celldatadir.c_str());
            cellgravesdatadir = strdup(name_cellgravedatadir.c_str());
            fooddatadir = strdup(name_fooddatadir.c_str());
            // this took a while to code :P
        } else {
            cerr << "Something went wrong reading the commandline arguments" << endl;
            return 1;
        }
    }
    return 0;
}

void Parameter::Read(const char *filename) {

    static bool ReadP = false;

    if (ReadP) {

        //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
        CleanUp();

    } else
        ReadP = true;

    FILE *fp = OpenReadFile(filename);


    T = fgetpar(fp, "T", T, true);
    target_area = igetpar(fp, "target_area", target_area, true);
    lambda = fgetpar(fp, "lambda", lambda, true);
    conn_diss = igetpar(fp, "conn_diss", conn_diss, true);
    border_energy = igetpar(fp, "border_energy", border_energy, true);
    neighbours = igetpar(fp, "neighbours", neighbours, true);
    min_area_for_life = igetpar(fp, "min_area_for_life", min_area_for_life, true);
    Jmed = fgetpar(fp, "Jmed", Jmed, true);
    Jalpha = fgetpar(fp, "Jalpha", Jalpha, true);
    n_chem = igetpar(fp, "n_chem", n_chem, true);
    if (n_chem) {
        diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
        decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
        secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
        dt = fgetpar(fp, "dt", dt, true);
        dx = fgetpar(fp, "dx", dx, true);
    }
    n_init_cells = igetpar(fp, "n_init_cells", n_init_cells, true);
    size_init_cells = igetpar(fp, "size_init_cells", size_init_cells, true);
    sizex = igetpar(fp, "sizex", sizex, true);
    sizey = igetpar(fp, "sizey", sizey, true);
    mcs = igetpar(fp, "mcs", mcs, true);
    rseed = ugetpar(fp, "rseed", rseed, true);
    save_movie_period = igetpar(fp, "save_movie_period", save_movie_period, true);
    graphics = bgetpar(fp, "graphics", graphics, true);
    save_movie = bgetpar(fp, "save_movie", save_movie, true);
    genomefile = sgetpar(fp, "genomefile", genomefile, true);
    nodivisions = bgetpar(fp, "nodivisions", nodivisions, true);
    nodeaths = bgetpar(fp, "nodeaths", nodeaths, true);
    nr_regnodes = igetpar(fp, "nr_regnodes", nr_regnodes, true);
    mu = fgetpar(fp, "mu", mu, true);
    mustd = fgetpar(fp, "mustd", mustd, true);
    divtime = igetpar(fp, "divtime", divtime, true);
    circle_thickness = igetpar(fp, "circle_thickness", circle_thickness, true);
    circle_segments = igetpar(fp, "circle_segments", circle_segments, true);
    circle_dist = igetpar(fp, "circle_dist", circle_dist, true);
    divdur = igetpar(fp, "divdur", divdur, true);
    gompertz_alpha = fgetpar(fp, "gompertz_alpha", gompertz_alpha, true);
    gompertz_beta = fgetpar(fp, "gompertz_beta", gompertz_beta, true);
    scatter_start = bgetpar(fp, "scatter_start", scatter_start, true);
    chemgrad = bgetpar(fp, "chemgrad", chemgrad, true);
    groupextinction = bgetpar(fp, "groupextinction", groupextinction, true);
    randgroups = bgetpar(fp, "randgroups", randgroups, true);
    dynamic_adh = bgetpar(fp, "dynamic_adh", dynamic_adh, true);
    chemcircles = bgetpar(fp, "chemcircles", chemcircles, true);
    miglines = bgetpar(fp, "miglines", miglines, true);
    existing_dirs = bgetpar(fp, "existing_dirs", existing_dirs, true);
    replace_dirs = bgetpar(fp, "replace_dirs", replace_dirs, true);
    moviedir = sgetpar(fp, "moviedir", moviedir, true);
    celldatafile = sgetpar(fp, "celldatafile", celldatafile, true);
    fooddatafile = sgetpar(fp, "fooddatafile", fooddatafile, true);
    save_data_period = igetpar(fp, "save_data_period", save_data_period, true);
    initial_food_amount = igetpar(fp, "initial_food_amount", initial_food_amount, true);
    metabperiod = igetpar(fp, "metabperiod", metabperiod, true);
    gradnoise = fgetpar(fp, "gradnoise", gradnoise, true);
    gradscale = fgetpar(fp, "gradscale", gradscale, true);
    foodpatches = igetpar(fp, "foodpatches", foodpatches, true);
    foodpatchperiod = igetpar(fp, "foodpatchperiod", foodpatchperiod, true);
    seasonduration = igetpar(fp, "seasonduration", seasonduration, true);
    seasonamplitude = fgetpar(fp, "seasonamplitude", seasonamplitude, true);
    foodpatcharea = igetpar(fp, "foodpatcharea", foodpatcharea, true);
    foodperspot = igetpar(fp, "foodperspot", foodperspot, true);
    maxfoodpatches = igetpar(fp, "maxfoodpatches", maxfoodpatches, true);
    foodstart = igetpar(fp, "foodstart", foodstart, true);
    eatperiod = igetpar(fp, "eatperiod", eatperiod, true);
    chancemediumcopied = fgetpar(fp, "chancemediumcopied", chancemediumcopied, true);
    colortablefile = sgetpar(fp, "colortablefile", colortablefile.c_str(), true);
    plots = sgetpar(fp, "plots", plots.c_str(), true);
    input_scales = sgetpar(fp, "input_scales", input_scales.c_str(), true);
    key_lock_weights = sgetpar(fp, "key_lock_weights", key_lock_weights.c_str(), true);
    circle_coverage = fgetpar(fp, "circle_coverage", circle_coverage, true);
    persduration = igetpar(fp, "persduration", persduration, true);
    startmu = fgetpar(fp, "startmu", startmu, true);
    init_chemmu = fgetpar(fp, "init_chemmu", init_chemmu, true);
    grn_update_period = igetpar(fp, "grn_update_period", grn_update_period, true);
    latticedir = sgetpar(fp, "latticedir", latticedir, true);
    celldatadir = sgetpar(fp, "celldatadir", celldatadir, true);
    cellgravesdatadir = sgetpar(fp, "cellgravesdatadir", cellgravesdatadir, true);
    fooddatadir = sgetpar(fp, "fooddatadir", fooddatadir, true);
    save_lattice_period = igetpar(fp, "save_lattice_period", save_lattice_period, true);
    evolreg = bgetpar(fp, "evolreg", evolreg, true);
    info_period = igetpar(fp, "info_period", info_period, true);
    init_cell_config = igetpar(fp, "init_cell_config", init_cell_config, true);
    cell_placement = igetpar(fp, "cell_placement", cell_placement, true);
}

// In the future the parser for the rules for key to J val tau,medium
// will be more developed, maybe even evolvable 8O
// int Parameter::SumLookupTableValue(int *lookup_table){
//   return -1;
// }
// int Parameter::MultiplyLookupTableValue(int *lookup_table){
//   return -1;
// }


const char *sbool(const bool &p) {

    const char *true_str = "true";
    const char *false_str = "false";
    if (p)
        return true_str;
    else
        return false_str;
}

void Parameter::Write(ostream &os) const {
    setlocale(LC_NUMERIC, "C");

    os << " T = " << T << endl;
    os << " target_area = " << target_area << endl;
    os << " lambda = " << lambda << endl;
    //if (Jtable)
    //  os << " Jtable = " << Jtable << endl;
    os << " conn_diss = " << conn_diss << endl;
    os << " border_energy = " << border_energy << endl;
    os << " neighbours = " << neighbours << endl;
    os << " min_area_for_life = " << min_area_for_life << endl;
    os << " Jmed = " << Jmed << endl;
    os << " Jalpha = " << Jalpha << endl;
    os << " n_chem = " << n_chem << endl;
    os << " diff_coeff = " << diff_coeff[0] << endl;
    os << " decay_rate = " << decay_rate[0] << endl;
    os << " secr_rate = " << secr_rate[0] << endl;
    os << " dt = " << dt << endl;
    os << " dx = " << dx << endl;
    os << " n_init_cells = " << n_init_cells << endl;
    os << " size_init_cells = " << size_init_cells << endl;
    os << " sizex = " << sizex << endl;
    os << " sizey = " << sizey << endl;
    os << " mcs = " << mcs << endl;
    os << " rseed = " << rseed << endl;
    os << " save_movie_period = " << save_movie_period << endl;
    os << " graphics = " << sbool(graphics) << endl;
    os << " save_movie = " << sbool(save_movie) << endl;
    if (genomefile) {
        os << " genomefile = " << genomefile << endl;
    }
    os << " nr_regnodes = " << nr_regnodes << endl;
    os << " mu = " << mu << endl;
    os << " mustd = " << mustd << endl;
    os << " divtime= " << divtime << endl;
    os << " circle_thickness= " << circle_thickness << endl;
    os << " circle_dist= " << circle_dist << endl;
    os << " circle_segments= " << circle_segments << endl;
    os << " divdur = " << divdur << endl;
    os << " gompertz_alpha = " << gompertz_alpha << endl;
    os << " gompertz_beta = " << gompertz_beta << endl;
    os << " initial_food_amount = " << initial_food_amount << endl;
    os << " metabperiod = " << metabperiod << endl;
    os << " gradnoise = " << gradnoise << endl;
    os << " gradscale = " << gradscale << endl;
    os << " foodpatches = " << foodpatches << endl;
    os << " foodpatchperiod = " << foodpatchperiod << endl;
    os << " seasonduration = " << seasonduration << endl;
    os << " seasonamplitude = " << seasonamplitude << endl;
    os << " foodpatcharea = " << foodpatcharea << endl;
    os << " foodperspot = " << foodperspot << endl;
    os << " maxfoodpatches = " << maxfoodpatches << endl;
    os << " foodstart = " << foodstart << endl;
    os << " eatperiod = " << eatperiod << endl;
    os << " nodivisions = " << nodivisions << endl;
    os << " nodeaths = " << nodeaths << endl;
    os << " chancemediumcopied = " << chancemediumcopied << endl;
    os << " celldatafile = " << celldatafile << endl;
    os << " fooddatafile = " << fooddatafile << endl;
    os << " save_data_period = " << save_data_period << endl;
    os << " colortablefile = " << colortablefile << endl;
    os << " plots = " << plots << endl;
    os << " key_lock_weights = " << key_lock_weights << endl;
    os << " input_scales = " << input_scales << endl;
    os << " circle_coverage = " << circle_coverage << endl;
    os << " persduration = " << persduration << endl;
    os << " startmu = " << startmu << endl;
    os << " grn_update_period = " << grn_update_period << endl;
    os << " latticedir = " << latticedir << endl;
    os << " celldatadir = " << celldatadir << endl;
    os << " cellgravedatadir = " << cellgravesdatadir << endl;
    os << " fooddatadir = " << fooddatadir << endl;
    os << " save_lattice_period = " << save_lattice_period << endl;
    if (moviedir)
        os << " moviedir = " << moviedir << endl;
    os << " evolreg = " << evolreg << endl;
    os << " info_period = " << info_period << endl;
    os << " init_cell_config = " << init_cell_config << endl;
    os << " cell_placement = " << cell_placement << endl;
}

ostream &operator<<(ostream &os, Parameter &p) {
    p.Write(os);
    return os;
}

Parameter par;
