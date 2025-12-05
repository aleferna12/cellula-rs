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


/// @file CompPolygon.h
#pragma once
#include "Edge2.h"
#include "common.h"

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif


/** \brief Connected component with boundary- and hole polygons
 *
 * The CompPolygon struct holds a connected component of triangles.
 * Thereby two triangles count as connected if they share a common
 * edge. It also holds an ordered vector of polygon edges and
 * ordered vectors of hole edges.
 *
 * @remark An edge is represented by a triangle and an opposite
 * index. Its orientation is always counterclockwise (CCW) around
 * the triangle. Thus outer polygon edges are counterclockwise
 * while hole polygons (if any) are clockwise (CW) oriented.
 */
struct CLASS_DECLSPEC CompPolygon
{
	CompPolygon();
	CompPolygon(const CompPolygon& other);
	CompPolygon& operator=(const CompPolygon& other);
	~CompPolygon();
	// Data
	std::vector<Triangle2*>* pvCC; ///< Connected component of triangles (connection along edges, not just vertices)
	std::vector<Edge2>* pvOuterPolygon; ///< Ordered outer polygon
	std::vector<std::vector<Edge2> >* pvvHolePolygons; ///< Ordered hole polygons
};


} // (namespace)
