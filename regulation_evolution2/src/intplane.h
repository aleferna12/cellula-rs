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

#ifndef _INTPL_HH_
#define _INTPL_HH_

#include <stdio.h>
#include <float.h>
#include <functional>
#include <climits>
#include "graph.h"
#include "random.h"
#include "foodpatch.h"

class IntPlane; //forward declaration

class CellularPotts;

class IntPlane {

    friend class Info;

public:

    /*!
     * \brief Constructor for PDE object.
     * \param sx: horizontal size of PDE planes
     * \param sy: vertical size of PDE planes
    */

    IntPlane(int sx, int sy, int fill = 0);

    // destructor must also be virtual
    virtual ~IntPlane();

    //! \brief Returns the horizontal size of the PDE planes.
    inline int SizeX() const {
        return sizex;
    }

    //! \brief Returns the vertical size of the PDE planes.
    inline int SizeY() const {
        return sizey;
    }

    double getDiagonal() const {
        return sqrt(sizex * sizex + sizey * sizey);
    }

    int getArea() const {
        return sizex * sizey;
    }

    /*! \brief Returns the value of grid point x,y of PDE plane "layer".

    Warning, no range checking done.

    \param x, y: grid point to probe.
    */
    inline int Sigma(const int x, const int y) const {
        return sigma[x * sizey + y];
    }

    /*! \brief Sets grid point x,y of PDE plane "layer" to value "value".

    \param x, y: grid point
    \param value: new contents

    */
    inline void setSigma(const int x, const int y, const int value) {
        sigma[x * sizey + y] = value;
    }

    int SetNextVal(int sig);

    tuple<int, int> getMinMax() const {
        int maxval = 0;
        int minval = INT_MAX;
        for (int i = 1; i < sizex - 1; ++i)
            for (int j = 1; j < sizey - 1; ++j) {
                minval = min(Sigma(i, j), minval);
                maxval = max(Sigma(i, j), maxval);
            }
        return make_tuple(minval, maxval);
    }

protected:
    // Protected member functions

    /*! \brief Used in makePlots. Takes a color and turns it into a grey value.

    \param val: Value from PDE plane.

    Implement this function in you main simulation code. See e.g. vessel.cpp.
    */
    //virtual int MapColour(double val);

    //! empty constructor (necessary for derivation)
    IntPlane();

private:
    int *sigma;
    int sizex;
    int sizey;
};

#endif
