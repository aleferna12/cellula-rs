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

/// @file Visualizer2.h
#pragma once
#include "Point2.h"
#include "Circle2.h"
#include "PShape.h"
#include "Segment2.h"
#include "Color.h"
#include "Label.h"
#include "Bbox2.h"
#include "Edge2.h"
#include "VertexPair2.h"


#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif
class ConstraintSegment2; // FWD
class ConstraintGraph2; // FWD
class VoroCell2; // FWD
struct Dat; // FWD
class SegmentChecker; // FWD

/** \brief Visualizer2 is a PDF- and Postscript writer.
 *
* \sa http://www.geom.at/example2-traversing/
* \image html visualizer.jpg "Example output of the Visualizer"
* \image latex visualizer.eps "Example output of the Visualizer" width=12cm
*
*
*/

class Visualizer2
{
public:
/** \brief Constructor
*
* @param filename_ must end with ".ps" or ".pdf".
*/
	CLASS_DECLSPEC
	explicit Visualizer2(const char* filename_);

	CLASS_DECLSPEC
	~Visualizer2();
	/** \brief Add a vector of Voronoi Cells to the visualization
	 *
	 * @param vC Vector of pointers to Voronoi cells.
	 * @param c Color for the Voronoi cells.  Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the Voronoi cells.
	*/
	CLASS_DECLSPEC
	void addObject(const std::vector<VoroCell2*>& vC,const Color& c);

	/** \brief Add polygonal shapes
	 *
	 * @param vPolygonalShapes Vector of polygonal shapes.
	 * @param c Color for the polygonal shapes. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the shapes.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<PShape>& vPolygonalShapes,const Color& c);

	/** \brief Add a polygonal shape
	 *
	 * @param polygonalShape A single polygonal shape.
	 * @param c Color for the polygonal shapes. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the shape.
	 */
	CLASS_DECLSPEC
	void addObject(const PShape& polygonalShape,const Color& c);

	/** \brief Add a Voronoi cell to the visualization
	 *
	 * @param pVoroCell Pointer to a Voronoi cell.
	 * @param c Color for the Voronoi cell. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the Voronoi cell.
	 */
	CLASS_DECLSPEC
	void addObject(VoroCell2* pVoroCell,const Color& c);

	/** \brief Add a Segment2 object to the visualization
	 *
	 * @param seg Segment2 object.
	 * @param c Color for the segment. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to draw marks at the endpoints.
	 *
	 */
	CLASS_DECLSPEC
	void addObject(const Segment2& seg,const Color& c);

	/** \brief Add a ConstraintGraph2 object to the visualization
	 *
	 * @param pCG Pointer to a ConstraintGraph2
	 * @param c Color for the ConstraintGraph2. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to draw marks at the endpoints.
	 *
	 */
	CLASS_DECLSPEC
	void addObject(ConstraintGraph2* pCG,const Color& c);


	/** \brief Add an Edge2 object to the visualization
	 *
	 * @param edge Edge2 object.
	 * @param c Color for the edge. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to draw marks at the endpoints.
	 */
	CLASS_DECLSPEC
	void addObject(const Edge2& edge,const Color& c);

	/** \brief Add a vector of Point2 objects to the visualization
	 *
	 * @param vPoints Vector of Point2 objects.
	 * @param c Color for the points.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Point2>& vPoints,const Color& c);

	/** \brief Add a vector of Point2 pointers to the visualization
	 *
	 * @param vPoints Vector of pointers to Point2 objects.
	 * @param c Color for the points.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Point2*>& vPoints,const Color& c);


	/** \brief Add a vector of Segment2 objects to the visualization
	 *
	 * @param vSegments Vector of Segment2 objects.
	 * @param c Color for the segments.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Segment2>& vSegments,const Color& c);

	/** \brief Add a vector of VertexPair2 objects as segments to the visualization
	 *
	 * @param vVertexPairs Vector of VertexPair2 objects.
	 * @param c Color for the segments.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<VertexPair2>& vVertexPairs,const Color& c);

	/** \brief Add a vector of ConstraintSegment2 pointers to the visualization
	 *
	 * @param vConstraintSegments Vector of pointers to ConstraintSegment2 objects.
	 * @param c Color for the constraint segments. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to draw marks at the endpoints.
	 */


	CLASS_DECLSPEC
	void addObject(const std::vector<ConstraintSegment2*>& vConstraintSegments,const Color& c);

	/** \brief Add a vector of Edge2 objects to the visualization
	 *
	 * @param vSegments Vector of Edge2 objects.
	 * @param c Color for the edges. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to draw marks at the endpoints.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Edge2>& vSegments,const Color& c);


	/** \brief Add a vector of Triangle2 objects to the visualization
	 *
	 * @param vT Vector of Triangle2 objects.
	 * @param c Color for the triangles. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the triangles
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Triangle2>& vT,const Color& c);

	/** \brief Add a Circle2 object to the visualization
	 *
	 * @param circ Circle2 object.
	 * @param c Color for the circle. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the circle
	 */
	CLASS_DECLSPEC
	void addObject(const Circle2& circ,const Color& c);

	/** \brief Add a Point2 object to the visualization
	 *
	 * @param pnt Point2 object.
	 * @param c Color for the point.
	 */	CLASS_DECLSPEC
	void addObject(const Point2& pnt,const Color& c);

	/** \brief Add a Triangle2 object to the visualization
	 *
	 * @param tri Triangle2 object.
	 * @param c Color for the triangle. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the triangle.
	 */
	CLASS_DECLSPEC
	void addObject(const Triangle2& tri,const Color& c);

	/** \brief Add a vector of Triangle2 pointers to the visualization
	 *
	 * @param vT Vector of pointers to Triangle2 objects.
	 * @param c Color for the triangles. Define the Color with bFill=true, i.e., Color(CRED,0.01,true), to fill the triangles.
	 */
	CLASS_DECLSPEC
	void addObject(const std::vector<Triangle2*>& vT,const Color& c);

	/** \brief Add a Label object to the visualization
	 *
	 * @param lab Label object.
	 * @param c Color for the label.
	 */
	CLASS_DECLSPEC
	void addObject(const Label& lab,const Color& c);

	/** \brief Add a header line to the visualization
	 *
	 * @param s Header line as a C-string.
	 */
	CLASS_DECLSPEC
	void addHeaderLine(const char* s);


	/** \brief Finish and write the output file
	 *
	 * @note This method \e must be called at the end when all the objects have been added.
	 */
	CLASS_DECLSPEC
	void writeFile();


/** @private
 */
	CLASS_DECLSPEC
	void setLimit(const Bbox2& bbx);

	/** \brief Compute the range
	 *
	 * @param bWithVoronoi specifies if the Voronoi cells shall be
	 * incorporated.
	 *
	 * @return a bounding box of currently contained objects.
	 */
	Bbox2 computeRange(bool bWithVoronoi);
protected:
	Visualizer2(const Visualizer2& );
	Visualizer2& operator=(const Visualizer2&);
	// Helpers
	Point2 scaledPoint(const Point2 &p);
	double scaledDouble(const double &d);
	void changeColor(float r,float g,float b,float linewidth,bool bFill);
	void changeColor(const Color& c);
	// Write Header/Footer
	void writeHeader(const char* title);
	void writeFooter();
	void writeHeaderLines();
	// Write Objects
	void writeLabel(Label l);
	void writeLine(const Point2& pSource,const Point2& pTarget);
	void writePolygon(PShape& pshape,bool bFill,double width);
	void writeTriangle(const Point2& p0_,const Point2& p1_,const Point2& p2_,bool bFill,double width);
	void writeTriangle(const Triangle2* pT,bool bFill_,double width);
	void writeVoroCell(VoroCell2* pVoroCell,bool bFill,double width);
	void writePoint(const Point2& p1_,float size);
	void writeMark(const Point2& p1_,float size);
	void writeCircle(const Point2& p1_,double radius,bool bFill);
	void periodicStroke();

	// DATA
	bool bTypePS;
	Dat* pDat;
	int updateCtr;
	bool bFill;
	Color lastColor;
	bool bFileClosed;
	Bbox2 bbox;
	std::ofstream outFile;
	std::vector<std::pair<Segment2,Color> > vSegments;
	std::vector<std::pair<Circle2,Color> > vCircles;
	std::vector<std::pair<Point2,Color> > vPoints;
	std::vector<std::pair<Triangle2,Color> > vTriangles;
	std::vector<std::pair<Label,Color> > vLabels;
	std::vector<std::pair<VoroCell2*,Color> > vVoroCells;
	std::vector<std::pair<PShape,Color> > vPolygons;
};


} // (namespace)
