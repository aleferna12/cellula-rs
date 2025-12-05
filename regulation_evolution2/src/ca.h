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

// mainpage.h contains no C++ code, it is for the main page of the
// documentation
#include "mainpage.h"

/*! Implementation of the Glazier & Graner cellular Potts model **/
#ifndef _CA_HH_
#define _CA_HH_

#include <vector>
#include <stdio.h>
#include <unordered_map>
#include <array>
#include "graph.h"
#include "pde.h"
#include "cell.h"
#include "boundingbox.h"

class Dish;

namespace std {
    //templates: typename is the variable type, which we can call T (can be int, char, vector, whatever)
    template<typename T, size_t N>
    struct hash<const array<T, N> > {
        typedef const array<T, N> argument_type;
        typedef size_t result_type;

        result_type operator()(const argument_type &a) const {
            hash<T> hasher;
            result_type h = 1;
            for (result_type i = 0; i < N; ++i) {
                h = h * hasher(a[i]);
            }
            return h;
        }
    };

    template<typename T, size_t N>
    struct hash<array<T, N> > {
        typedef array<T, N> argument_type;
        typedef size_t result_type;

        result_type operator()(const argument_type &a) const {
            hash<T> hasher;
            result_type h = 1;
            for (result_type i = 0; i < N; ++i) {
                h = h * hasher(a[i]);
            }
            return h;
        }
    };
}

//template <class T>
class Edges {
public:
    // this defines the operator [] for the this unordered_map (in the const and not const case <---? why the latter though)
    const std::unordered_map<std::array<int, 4>, int>::iterator &
    operator[](unsigned index) const { return edgevector[index]; }

    std::unordered_map<std::array<int, 4>, int>::iterator &operator[](unsigned index) { return edgevector[index]; }

    //Edges(){};
    //~Edges(){};

    bool insert(std::array<int, 4> edge) {
        int n = int(edgevector.size());
        auto insertion = edgemap.insert({edge, n});
        //unordered_map.insert()returns a pair object. first el is iterator to new element and second is bool indicating success of insert
        if (insertion.second) {
            edgevector.push_back(insertion.first); //store iterator pointing to new element in edgemap.
        }
        return insertion.second;
    }

    bool erase(std::array<int, 4> edge) {
        int index = edgemap.at(edge);
        // std::unordered_map<std::array<int,4>,int>::iterator last_it = edgevector
        std::swap(edgevector[index],
                  edgevector.back()); //swap the item to be removed with the item at the back, then pop out of the vector
        edgevector[index]->second = index; //item that used to be at the back gets the other's index
        edgevector.pop_back();
        edgemap.erase(edge);
        return true;
    }

    unsigned size_map() const {
        return edgemap.size();
    }

    unsigned size_vector() const {
        return edgevector.size();
    }

private:
    std::unordered_map<std::array<int, 4>, int> edgemap; //the edges are the keys, the value is length of edgevector upon insertion.
    std::vector<std::unordered_map<std::array<int, 4>, int>::iterator> edgevector; //iterators to elements in edgemap, ordered by time of insertion
};

class Dir {

    /* To store a celldirection matrix */
    friend class CellularPotts;

public:

    Dir() {
        meanx = -1;
        meany = -1.;
        aa1 = 0.;
        aa2 = 0.;
        bb1 = 0.;
        bb2 = 0.;
        lb1 = 0.;
        lb2 = 0.;
    }

    double meanx, meany;
    double aa1, aa2;
    double bb1, bb2;
    double lb1, lb2;
};


class CellularPotts {
public:
    //! \brief Constructs a CA field. This should be done in "Dish".
    explicit CellularPotts(std::vector<Cell> *cells, int sizex = 200,
                           int sizey = 200);

    // empty constructor
    // (necessary for derivation)
    CellularPotts();

    // Keyword virtual means, that derived classed (cppvmCellularPotts) can override
    // this function and carry out the memory allocation in their preferred way
    // Every time AllocateSigma is called in the base class methods
    // the function belonging the actual type will be called
    virtual void AllocateSigma(int sx, int sy);

    // destructor must also be virtual
    virtual ~CellularPotts();

    void InitializeEdgeList(bool init_single_cell_center); //Set the initial edgelist which are eligible to change

    //! Return the total area occupied by the cells
    inline int Mass() {
        int mass = 0;
        for (int i = 0; i < sizex * sizey; i++) {
            if (sigma[0][i] > 0) mass++;
        }
        return mass;
    }

    // Tells the cell cell_sigma to divide, given its boundingbox
    // TODO: Put in cell.h
    int DivideCell(int cell_sigma, BoundingBox box);

    /*! Implements the core CPM algorithm. Carries out one MCS.
      \return Total energy change during MCS.
    */
    int AmoebaeMove2(PDE *PDEfield);


    // int BlobCounting(void); (not implemented?)

    // Used internally to assign a Cell object to each "blob"
    // "blobs" should already have different indices.
    //
    // (i.e. to start from a binary image you'd probably first want
    // to apply a blob counting algorithm
    void ConstructInitCells(Dish &beast);

    //! Returns the number of completed Monte Carlo steps.
    inline int Time() const {
        return thetime;
    }


    // not currently used? In Critter implementation (see Hogeweg
    // 2000) this was used to have cells divide at double their original area.
    inline int ZygoteArea() const {
        return zygote_area;
    }

    //! \brief Return the horizontal size of the CA plane.
    inline int SizeX() const {
        return sizex;
    }

    //! \brief Return the vertical size of the CA plane.
    inline int SizeY() const {
        return sizey;
    }

    /*! \brief Return the value of lattice site (x,y).

    i.e. This will return the index of the cell which occupies site (x,y). */
    inline int Sigma(const int x, const int y) const {
        return sigma[x][y];
    }

    //this function facilitates creation from a backup file
    int SetNextSigma(int sig);

    /*! In this method the principal axes of the cells are computed using
     the method described in "Biometry", box 15.5
     \return a pointer to a "new[]"ed array containing the directions.
     The memory has to be freed afterwards using the delete[] operator
    */
    Dir *FindCellDirections3() const;

    int GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y);

    int PlaceCellsRandomly(int n_cells, int cellsize);

    int PlaceCellsOrderly(int n_cells, int cellsize);

    //! \brief Adds a new Cell and returns a reference to it.
    inline Cell &AddCell(Dish &beast) {
        cell->push_back(Cell(beast));
        return cell->back();
    }

    void RemoveCell(Cell *cellid, int min_area, int meanx, int meany);

    int FancyloopX(int loopdepth, int meanx, int meany, int thissig, bool above);

    int FancyloopY(int loopdepth, int meanx, int meany, int thissig, bool left);

    inline Cell &getCell(int c) {
        return (*cell)[c];
    }


    //This used to be private but I want to use it elsewhere
    static const int nx[25], ny[25];
    static const int nbh_level[5];
    int n_nb;

    void MeasureCellSizes();

    void killCell(int c_sigma);

    int getKLStrength(int key_dec, int lock_dec) const {
        return KL_strengths[key_dec * max_key_lock_dec + lock_dec];
    }

    //! Creates a table in the heap containing the energy difference values for every possible combination of
    //! keys and locks
    static float * allocateKLStrengths(vector<float> &rule, int max_kl_dec, bool dynamic_adh);

    double energyDifference(int sigma1, int sigma2) const;

    double adhesionEnergy(int sigma1, int sigma2) const;

    // Not used in the simulations but useful to output data
    double calculateGamma(int sigma1, int sigma2) const;

    Edges &getEdgesVector() {
        return edgeSetVector;
    }

private:
    int DeltaH(int x, int y, int xp, int yp, PDE *PDEfield = nullptr);

    int DeltaHWithMedium(int x, int y, PDE *PDEfield = nullptr);

    void ConvertSpin(int x, int y, int xp, int yp);

    void ConvertSpinToMedium(int x, int y);

    static int CopyvProb(int DH, double stiff);

    void MeasureCellSize(Cell &c);

    static void CopyProb(double T);

    bool ConnectivityPreservedP(int x, int y);

    // little debugging function to print the site and its neighbourhood
    inline void PrintSite(int x, int y) {
        std::cout << "--------\n";
        std::cout << "[" << sigma[x - 1][y - 1] << " " << sigma[x][y - 1] << " " << sigma[x + 1][y - 1] << "]\n";
        std::cout << "[" << sigma[x - 1][y] << " " << sigma[x][y] << " " << sigma[x + 1][y] << "]\n";
        std::cout << "[" << sigma[x - 1][y + 1] << " " << sigma[x][y + 1] << " " << sigma[x + 1][y + 1] << "]\n";
    }


protected:
    void BaseInitialisation(std::vector<Cell> *cell);

    int **sigma;
    int sizex;
    int sizey;

private:
    Edges edgeSetVector;
    bool frozen;
    std::vector<Cell> *cell{};
    int zygote_area;
    int thetime;
    float *KL_strengths;

    int max_key_lock_dec;

    static float calculateKLStrength(unsigned long key_dec, unsigned long lock_dec, const vector<float> &rule) ;
};


#endif
