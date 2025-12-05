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

//This version of random.h of Tissue Simulation Toolkit
// uses c++11 <random> module, a giant step into modernity :P

#ifndef _RND_HH_
#define _RND_HH_

#include <random> //--- FOR THIS YOU NEED c++11, enable somehow with -std=c++11 flag
#include "parameter.h"

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// Declare engine

//extern std::mt19937_64 my_rng;
extern std::mt19937_64 my_rng;

//Declare distributions:
extern std::uniform_real_distribution<double> my_unif_real_dist;
// extern std::uniform_int_distribution<long> my_unif_int_dist; //not using this
extern std::binomial_distribution<int> my_binomial_dist;

unsigned int Seed(unsigned int seed);

double RANDOM();

double RANDNORMAL(double mean, double stdev);

long RandomNumber(long max);

void AskSeed();

int Randomize(void);

int BinomialDeviate(int N, double p); //returns a random number from binomial distr.
// with parameters N,p



#endif
