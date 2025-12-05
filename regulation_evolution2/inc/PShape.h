// Copyright (C) Geom Software e.U, Bernhard Kornberger, Graz/Austria
//
// This file is part of the Fade2D library. The student license is free
// of charge and covers personal non-commercial research. Licensees
// holding a commercial license may use this file in accordance with
// the Commercial License Agreement.
//
// This software is provided AS IS with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are
// not clear to you.
//
// Author: Bernhard Kornberger, bkorn (at) geom.at
// http://www.geom.at


/// @file PShape.h
#pragma once
#include "Point2.h"
#include "common.h"

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

/** \brief Polygonal Shape for Visualization
* \see Visualizer2
*/

class CLASS_DECLSPEC PShape
{
public:
	/** \brief Constructor
	*/
	explicit PShape(std::vector<Point2>& vP_);
	PShape(const PShape& other);
	PShape& operator=(const PShape& other);
	~PShape();
	friend std::ostream &operator<<(std::ostream &stream, PShape b);
protected:
	friend class Visualizer2;
	std::vector<Point2>* pVP;
};

} // (namespace)
