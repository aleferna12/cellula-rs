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

/// @file VtkWriter.h
#pragma once
#include "common.h"
#include "Point2.h"
#include "Segment2.h"
#include "Edge2.h"
#include "Triangle2.h"

#if GEOM_PSEUDO3D == GEOM_TRUE
namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D == GEOM_FALSE
namespace GEOM_FADE2D {
#else
#error GEOM_PSEUDO3D is not defined
#endif

struct Dat3; // FWD

/// \enum VtkColor
/// \brief Enumeration of colors used by the VTKWriter
enum VtkColor
{
	VTK_RED,
	VTK_GREEN,
	VTK_BLUE,
	VTK_YELLOW,
	VTK_CYAN,
	VTK_MAGENTA,
	VTK_BLACK,
	VTK_WHITE,
	VTK_GRAY,
	VTK_ORANGE,
	VTK_PINK,
	VTK_LIGHTBLUE,
	VTK_DARKBLUE,
	VTK_LIGHTGREEN,
	VTK_BROWN,
	VTK_PURPLE,
	VTK_TRANSPARENT ///< Meaning: Don't draw
};

/// \class VtkWriter
/// \brief A writer for the VTK file format
class CLASS_DECLSPEC VtkWriter
{
public:
	/// \brief Constructor
	/// \param name Name of the VTK file to write
	explicit VtkWriter(const char* name);

	/// \brief Destructor
	~VtkWriter();

	/// \brief Get the next color in the sequence
	/// \return The next VtkColor
	VtkColor getNextColor();

	/// \brief Add a point with a specified color
	/// \param p The point to add
	/// \param color The color of the point
	void addPoint(const Point2& p, VtkColor color);

	/// \brief Add multiple points with a specified color
	/// \param vPoints Vector of points to add
	/// \param color The color of the points
	void addPoints(const std::vector<Point2*>& vPoints, VtkColor color);
	void addPoints(const std::vector<Point2>& vPoints, VtkColor color);

	/// \brief Add a segment with a specified color
	/// \param segment The segment to add
	/// \param color The color of the segment
	/// \param bWithEndPoints Whether to include end points
	void addSegment(const Segment2& segment, VtkColor color, bool bWithEndPoints = false);

	/// \brief Add a segment defined by source and target points
	/// \param source, target Endpoints of the segment
	/// \param color The color of the segment
	/// \param bWithEndPoints Whether to include end points
	void addSegment(const Point2& source, const Point2& target, VtkColor color, bool bWithEndPoints = false);

	/// \brief Add multiple segments with a specified color
	/// \param vSegments Vector of segments to add
	/// \param color The color of the segments
	/// \param bWithEndPoints Whether to include end points
	void addSegments(const std::vector<Segment2>& vSegments, VtkColor color, bool bWithEndPoints = false);

	/// \brief Add multiple Edge2 objects with a specified color
	/// \param vEdges Vector of edges to add
	/// \param color The color of the edges
	/// \param bWithEndPoints Whether to include end points
	void addSegments(const std::vector<Edge2>& vEdges, VtkColor color, bool bWithEndPoints = false);

	/// \brief Add multiple segments with a specified color
	/// \param vSegmentEndPoints specifies n segments by 2*n endpoints
	/// \param color The color of the segments
	/// \param bWithEndPoints Whether to include end points
	void addSegments(const std::vector<Point2>& vSegmentEndPoints, VtkColor color, bool bWithEndPoints = false);

	/// \brief Add multiple triangles with a specified color
	/// \param vT Vector of triangle pointers to add
	/// \param color The color of the triangles
	void addTriangles(const std::vector<Triangle2*>& vT, VtkColor color);

	/// \brief Add multiple triangles defined by their corner points
	/// \param vTriangleCorners Vector of points defining the triangles
	/// \param color The color of the triangles
	void addTriangles(const std::vector<Point2>& vTriangleCorners, VtkColor color);

	/// \brief Add a triangle with a specified color
	/// \param t The triangle to add
	/// \param color The color of the triangle
	void addTriangle(const Triangle2& t, VtkColor color);

	/// \brief Add a triangle defined by three points
	/// \param p0, p1, p2 Points defining the triangle
	/// \param color The color of the triangle
	void addTriangle(const Point2& p0, const Point2& p1, const Point2& p2, VtkColor color);

	/// \brief Write the VTK file
	void writeFile();

private:
	/// \brief Copy constructor (deleted)
	VtkWriter(const VtkWriter&);

	/// \brief Assignment operator (deleted)
	VtkWriter& operator=(const VtkWriter&);

	Dat3* pDat; ///< Pointer to the data structure
};

} // namespace GEOM_FADE25D or GEOM_FADE2D

