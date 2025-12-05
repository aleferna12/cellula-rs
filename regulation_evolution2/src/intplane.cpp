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
#include <math.h>
#include <cstdlib>
#include <functional> // for bind, see InitIncreaseVal()
#include "crash.h"
#include "parameter.h"
#include "ca.h"
#include "intplane.h"
#include "misc.h"

using std::placeholders::_1;

extern Parameter par;

/** PRIVATE **/

IntPlane::IntPlane(const int sx, const int sy, int fill) {
    sizex = sx;
    sizey = sy;

    sigma = new int[sizex * sizey]{};
    if (fill != 0) {
        fill_n(sigma, sizex * sizey, fill);
    }
}


IntPlane::IntPlane() {
    sigma = nullptr;
    sizex = 0;
    sizey = 0;
}

// destructor (virtual)
IntPlane::~IntPlane() {
    delete[] sigma;
}

//copy of this function in ca.cpp
int IntPlane::SetNextVal(int pos) {
    //the plane has a 1 px boundary on all size, therefore we place the pixels
    //within that boundary
    static int xcount = 1, ycount = 1;

    if (xcount >= sizex - 1 || ycount >= sizey - 1) {
        return 1;
    }

    setSigma(xcount, ycount, pos);
    ycount++;
    if (ycount == sizey - 1) {
        ycount = 1;
        xcount++;
    }
    return 0;
}


