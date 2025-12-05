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

/// @file Zone2.h

#pragma once

#include "common.h"
#include "freeFunctions.h"
#include "FadeExport.h"
#include "Bbox2.h"
#include "Edge2.h"
#include "Segment2.h"
#include "UserPredicates.h"
#include "MsgBase.h"
#include "VtkWriter.h"
#include "CompPolygon.h"
#include "PolygonClipper.h"

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

class ZoneShooter; // FWD
class ConstraintSegment2; // FWD

/** \enum OptimizationMode
*
* Enumerates the possible modes for Valley/Ridge optimization
* through Zone2::slopeValleyRidgeOptimization().
*
*/

enum OptimizationMode
{
	OPTMODE_STANDARD, ///< Fastest optimization mode
	OPTMODE_BETTER, ///< Considerably better quality and still fast
	OPTMODE_BEST ///< Best quality but quite time consuming
};

class Progress; // FWD
class Dt2; // Fwd
class ConstraintGraph2; // Fwd
class Triangle2; // Fwd
class Point2; // Fwd
class Visualizer2; // Fwd
class Visualizer3; // Fwd

//
///** \brief Connected component with boundary- and hole polygons
// *
// * The CompPolygon struct holds a connected component of triangles.
// * Thereby two triangles count as connected if they share a common
// * edge. It also holds an ordered vector of polygon edges and
// * ordered vectors of hole edges.
// *
// * @remark An edge is represented by a triangle and an opposite
// * index. Its orientation is always counterclockwise (CCW) around
// * the triangle. Thus outer polygon edges are counterclockwise
// * while hole polygons (if any) are clockwise (CW) oriented.
// */
//struct CLASS_DECLSPEC CompPolygon
//{
//	CompPolygon();
//	CompPolygon(const CompPolygon& other);
//	CompPolygon& operator=(const CompPolygon& other);
//	~CompPolygon();
//	// Data
//	std::vector<Triangle2*>* pvCC; ///< Connected component of triangles (connection along edges, not just vertices)
//	std::vector<Edge2>* pvOuterPolygon; ///< Ordered outer polygon
//	std::vector<std::vector<Edge2> >* pvvHolePolygons; ///< Ordered hole polygons
//};

/** \brief Zone2 is a certain defined area of a triangulation
*
* \sa http://www.geom.at/example4-zones-defined-areas-in-triangulations/
* \sa http://www.geom.at/boolean-operations/
* \sa \ref createZone in the Fade2D class
*/
class CLASS_DECLSPEC Zone2
{
public:
	/// @private
	Zone2(Dt2* pDt_,ZoneLocation zoneLoc_);
	/// @private
	Zone2(Dt2* pDt_,ZoneLocation zoneLoc_,ConstraintGraph2* pConstraintGraph_);
	/// @private
	Zone2(Dt2* pDt_,ZoneLocation zoneLoc_,const std::vector<ConstraintGraph2*>& vConstraintGraphs_);
	/// @private
	Zone2(Dt2* pDt_,const std::vector<ConstraintGraph2*>& vConstraintGraphs_,ZoneLocation zoneLoc_,std::vector<Point2>& vStartPoints);
	/// @private
	~Zone2();

	/**
	* @brief Checks if the given point lies on the zone or its boundary.
	*
	* @param p The point to check.
	* @return true if the point is on the zone or its boundary, false otherwise.
	*
	* The first call to Zone2::hasOn(), Zone2::hasOnBoundary(), Zone2::getNearbyBoundaryConstraint(),
	* Zone2::shiftToZone() or Zone2::locate() initializes a search structure. For
	* zones not of type ZL_INSIDE, ZL_BOUNDED, or ZL_GLOBAL, constraint edges are inserted
	* around the zone.
	*/
	bool hasOn(const Point2& p);

	/**
	* @brief Checks if the given point lies on the boundary of the zone
	*
	* @param p The point to check.
	* @return true if the point is on the boundary, false otherwise.
	*
	* The first call to Zone2::hasOn(), Zone2::hasOnBoundary(), Zone2::getNearbyBoundaryConstraint(),
	* Zone2::shiftToZone() or Zone2::locate() initializes a search structure. For
	* zones not of type ZL_INSIDE, ZL_BOUNDED, or ZL_GLOBAL, constraint edges are inserted
	* around the zone.
	*/
	bool hasOnBoundary(const Point2& p);

	/**
	* @brief Locates the triangle containing the given point within the zone.
	*
	* @param p The point to locate.
	* @return A pointer to the Triangle2 if the point is inside the zone, otherwise NULL.
	*
	* The first call to Zone2::hasOn(), Zone2::hasOnBoundary(), Zone2::getNearbyBoundaryConstraint(),
	* Zone2::shiftToZone() or Zone2::locate() initializes a search structure. For
	* zones not of type ZL_INSIDE, ZL_BOUNDED, or ZL_GLOBAL, constraint edges are inserted
	* around the zone.
	*/
	Triangle2* locate(const Point2& p);

	/**
	* @brief Locates the nearest boundary ConstraintSegment2 of the zone within a specified distance.
	*
	* @param p The point from which to search.
	* @param tolerance The maximum allowed 2D distance between the point and the ConstraintSegment2.
	*
	* @return A pointer to the closest alive ConstraintSegment2 within the tolerance distance,
	* where p has an orthogonal projection, or NULL if no such segment is found.
	*
	* The first call to Zone2::hasOn(), Zone2::hasOnBoundary(), Zone2::getNearbyBoundaryConstraint(),
	* Zone2::shiftToZone() or Zone2::locate() initializes a search structure. For
	* zones not of type ZL_INSIDE, ZL_BOUNDED, or ZL_GLOBAL, constraint edges are inserted
	* around the zone.
	*/
	ConstraintSegment2* getNearbyBoundaryConstraint(Point2& p, double tolerance);
	/**
	* @brief Finds a point close to the input point that lies inside the zone.
	*
	* @param from A point outside the zone.
	* @param tolerance The maximum allowed 2D distance to find a point within the zone.
	* @param result A reference to store the point found within the zone, close to `from`.
	* @return true if a point within the zone is successfully found or created; false otherwise.
	*
	* The first call to Zone2::hasOn(), Zone2::hasOnBoundary(), Zone2::getNearbyBoundaryConstraint(),
	* Zone2::shiftToZone() or Zone2::locate() initializes a search structure. For
	* zones not of type ZL_INSIDE, ZL_BOUNDED, or ZL_GLOBAL, constraint edges are inserted
	* around the zone.
	*/
	bool shiftToZone(const Point2& from, double tolerance, Point2& result);


	/** \brief Save the zone
	 *
	 * This command saves the present Zone2 to a binary file. Any
	 * constraint edges and custom indices in the domain are retained.
	 *
	 * @param [in] filename is the output filename
	 * @return whether the operation was successful

	 * @note A Delaunay triangulation is convex without holes but this
	 * may not hold for the zone to be saved. Thus extra triangles may
	 * be saved to fill concavities. These extra-triangles will belong
	 * to the Fade_2D instance but not to the Zone2 object when reloaded.
	 *
	 * \sa save(std::ostream& stream). Use the similar command
	 * Fade_2D::saveZones(const char* file, std::vector<Zone2*>& vZones)
	 * to store more than just one zone. Use Fade_2D::saveTriangulation()
	 * to store all triangles of the triangulation plus any specified zones.
	 * Use Fade_2D::load() to reload the data from such files.
	*/
	bool save(const char* filename);

	/** \brief Save the zone
	 *
	 * This command saves the present Zone2 to an ostream. Any
	 * constraint edges and custom indices in the domain are retained.
	 *
	 * @param stream is the output stream
	 * @return whether the operation was successful

	 * @note A Delaunay triangulation is convex without holes but this
	 * may not hold for the zone to be saved. Thus extra triangles may
	 * be saved to fill concavities. These extra-triangles will belong
	 * to the Fade_2D instance but not to the Zone2 object when reloaded.
	 *
	 * \sa Use the similar command Fade_2D::saveZones(const char* file, std::vector<Zone2*>& vZones)
	 * to store more than just one zone. Use Fade_2D::saveTriangulation()
	 * to store all triangles of the triangulation plus any specified zones.
	 * Use Fade_2D::load() to reload the data from such files.
	*/
	bool save(std::ostream& stream);


/** \brief Export triangles from a zone
 *
 * @param fadeExport is a struct that will hold the requested triangulation data
 * @param bWithCustomIndices determines whether the custom indices of the points are also stored
 *
 */
	void exportZone(FadeExport& fadeExport,bool bWithCustomIndices) const;

/** \brief Get the boundary of an offset shape.
*
* This function computes a zone offset by a specified distance. The result is returned as
* a vector of oriented Segment2's.
*
* \param offset The positive or negative offset distance from the present zone.
* \param [out] vOffsetBoundary A vector that will contain the output segments, in arbitrary order, oriented counterclockwise around the resulting shape.
* \param mergeAngleDeg If the angle between two offset points of an original vertex is less than this
*                      value (in degrees), the points will be merged. Default: 10.0, valid range: greater than 0 and up to 135.
* \param angleStepDeg Specifies the angle interval (in degrees) at which circular arcs are sampled using line segments.
*                     Default: 20.0, valid range: greater than 0 and up to 135.
*/

	void getOffsetBoundary(double offset,std::vector<Segment2>& vOffsetBoundary,double mergeAngleDeg=10.0,double angleStepDeg=20.0);

	/** \brief Register a message receiver
	 *
	 * @param msgType is the type of message the subscriber shall receive, e.g. MSG_PROGRESS or MSG_WARNING
	 * @param pMsg is a pointer to a custom class derived from MsgBase
	*/
	void subscribe(MsgType msgType,MsgBase* pMsg);
	/** \brief Unregister a message receiver
	 *
	 * @param msgType is the type of message the subscriber shall not receive anymore
	 * @param pMsg is a pointer to a custom class derived from MsgBase
	*/
	void unsubscribe(MsgType msgType,MsgBase* pMsg);
	/** \brief Get the zone location
	*
	* \returns ZL_INSIDE if the zone applies to the triangles inside one or more ConstraintGraph2 objects@n
	* ZL_OUTSIDE if the zone applies to the outside triangles@n
	* ZL_GLOBAL if the zone applies (dynamically) to all triangles@n
	* ZL_RESULT if the zone is the result of a set operation@n
	* ZL_GROW if the zone is specified by a set of constraint graphs and an inner point@n
	* \image html in_and_outside_zone.jpg "An ouside zone and in inside zone"
	* \image latex in_and_outside_zone.eps "An ouside zone and in inside zone" width=12cm
	*/
	ZoneLocation getZoneLocation() const;

	/** \brief Convert a zone to a bounded zone
	*
	* \anchor convertToBoundedZone
	* The mesh generation algorithms refine() and refineAdvanced() require
	* a zone object that is bounded by constraint segments. This is always
	* the case for zones with zoneLocation ZL_INSIDE but other types of
	* zones may be unbounded. For convenience this method is provided to
	* create a bounded zone from a possibly unbounded one.
	*
	* @return a pointer to a new Zone2 object with zoneLocation ZL_RESULT_BOUNDED
	* or a pointer to the present zone if this->getZoneLocation() is ZL_INSIDE.
	*/
	Zone2* convertToBoundedZone();

	/** \brief Postscript- and PDF-visualization
	*
	* @param filename is the name of the output file.
	* @param bShowFull specifies if only the zone or the full triangulation shall be drawn
	* @param bWithConstraints specifies if constraint edges shall be drawn
	*
	*/
	void show(const char* filename,bool bShowFull,bool bWithConstraints) const;

	/** \brief Postscript- and PDF-visualization
	*
	* @param pVisualizer is a pointer to an existing Visualizer2 object drawing a .ps or .pdf file
	* @note You must call pVisualizer->writeFile() before program end
	* @param bShowFull specifies if only the zone or the full triangulation shall be drawn
	* @param bWithConstraints specifies if constraint edges shall be drawn
	*/
	void show(Visualizer2* pVisualizer,bool bShowFull,bool bWithConstraints) const;



	/** \brief VTK visualization
	*
	* @param filename The name of the output file.
	* @param zoneColor The color for the zone's triangles
	* @param nonZoneColor The color of the non-zone-triangles. Use VTK_TRANSPARENT to prevent them from being drawn.
	* @param constraintColor The color of the constraint edges. Use VTK_TRANSPARENT to prevent them from being drawn.
	*/
	void showVtk(const char* filename,VtkColor zoneColor,VtkColor nonZoneColor=VTK_TRANSPARENT,VtkColor constraintColor=VTK_TRANSPARENT) const;


	/** \brief VTK visualization
	*
	* @param pVtk A VtkWriter object that may already contain other geometric objects
	* @param zoneColor The color for the zone's triangles
	* @param nonZoneColor The color of the non-zone-triangles. Use VTK_TRANSPARENT to prevent them from being drawn.
	* @param constraintColor The color of the constraint edges. Use VTK_TRANSPARENT to prevent them from being drawn.
	*
	* @note pVtk->writeFile() finally writes the .vtk file.
	*/
	void showVtk(VtkWriter* pVtk,VtkColor zoneColor,VtkColor nonZoneColor=VTK_TRANSPARENT,VtkColor constraintColor=VTK_TRANSPARENT) const;



#if GEOM_PSEUDO3D==GEOM_TRUE
	/** \brief Geomview visualization
	*
	* @param filename is the name of the output file.
	* @param color is a string ("red green blue alpha"), e.g., "1.0 0.0 0.0 1.0"*
	*
	*/
	void showGeomview(const char* filename,const char* color) const;
	/** \brief Geomview visualization
	*
	* @param pVis points to a Visualizer3 object
	* @param color is a string ("red green blue alpha"), e.g., "1.0 0.0 0.0 1.0"*
	*
	*/
	void showGeomview(Visualizer3* pVis,const char* color) const;

	/// @private
	void analyzeAngles(const char* name="");

	/**
	 * \brief Smoothing
	 *
	 * **Deprecated method:** This method is deprecated and retained for backwards
	 * compatibility. It is recommended to use the new method @ref smoothing2()
	 * instead for better results.
	 *
	 * Weighted laplacian smoothing for the z-coordinate of all zone
	 * vertices. The x and y coordinates can also be optimized but
	 * only for vertices not belonging to the boundary of the zone
	 * or to constraint edges. This method is very fast but does
	 * nevertheless support the progress bar mechanism.
	 *
	 * @param numIterations is the number of smoothing passes.
	 * @param bWithXY specifies if the x and y coordinates are also adapted.
	 *
	*/
	void smoothing(int numIterations=2,bool bWithXY=true);

	/** \brief Smooths the vertices of the zone.
	 *
	 * This function performs a weighted Laplacian smoothing on the vertices of
	 * the Zone2. We distinguish between two types of vertices:
	 *  - **Static vertices**: Belong to constraint-edges or border-edges of the Zone2.
	 *  - **Dynamic vertices**: All other vertices.
	 *
	 * @param numIterations Number of smoothing iterations to perform. Higher
	 * values yield smoother results; 2-3 passes are recommended.
	 * @param bWithXY Boolean flag to determine if the x- and y-coordinates of dynamic
	 *                vertices should also be adjusted.
	 * @param bConstraintsWithZ Boolean flag to control whether the z-coordinates of
	 *                          static vertices are also updated by the smoothing process.
	 *
	 * @note
	 * - The x- and y-coordinates of static vertices are never changed.
	 * - The z-coordinates of dynamic vertices can always change, regardless of the
	 *   value of bConstraintsWithZ.
	 *
	 */
	void smoothing2(int numIterations,bool bWithXY,bool bWithConstraintZ);

	/** \brief Optimize Slopes, Valleys and Ridges
	 *
	 * A pure Delaunay triangulation takes only the x and y coordinates
	 * into account. However, for terrain scans, it is important to
	 * consider the z coordinate as well, otherwise ridges, valleys and
	 * rivers will look unnatural. This method leaves the points
	 * constant, but uses edge flips to change the connectivity, making
	 * the surface smoother overall.
	 *
	 * @param om is the optimization mode: OPTMODE_NORMAL is the fastest.
	 * OPTMODE_BETTER provides significantly better results while still
	 * taking a moderate amount of time. OPTMODE_BEST delivers the best
	 * results, but also has a significantly higher time requirement.
	 * This method supports the progress-bar mechanism.
	 *
	 * @note Flipping edges makes the triangulation non-delaunay, i.e.
	 * the empty-circle-property is then no longer given. Improving the
	 * smoothness of the surface by edge flips also means degrading the
	 * interior angles of the triangles (to a certain degree).
	 */
	void slopeValleyRidgeOptimization(OptimizationMode om=OPTMODE_BETTER);

	/// @private
	/*
	 * This function is deprecated but kept for backwards compatibility.
	 * Better use slopeValleyRidgeOptimization() (see above)
	 *
	 * Optimize Valleys and Ridges
	 *
	 * A Delaunay triangulation is not unique when when 2 or more triangles
	 * share a common circumcircle. As a consequence the four corners of
	 * a rectangle can be triangulated in two different ways: Either the
	 * diagonal proceeds from the lower left to the upper right corner
	 * or it connects the other two corners. Both solutions are valid and
	 * an arbitrary one is applied when points are triangulated. To improve
	 * the repeatability and for reasons of visual appearance this method
	 * unifies such diagonals such that they point from the lower left to
	 * the upper right corner (or in horizontal direction).\n
	 *
	 * Moreover a Delaunay triangulation does not take the z-value into
	 * account and thus valleys and ridges may be disturbed. The present
	 * method flips diagonals such that they point from the lower left to
	 * the upper right corner of a quad. And if the 2.5D lengths of the
	 * diagonals are significantly different, then the shorter one is
	 * applied.
	 *
	 * @param tolerance2D is 0 when only exact cases of more than 3 points
	 * on a common circumcircle shall be changed. But in practice input
	 * data can be disturbed by noise and tiny rounding errors such that
	 * grid points are not exactly on a grid. The numeric error is computed
	 * as \f$error=\frac{abs(diagonalA-diagonalB)}{max(diagonalA,diagonalB)}\f$.
	 * and \p tolerance2D is an upper threshold to allow modification despite
	 * such tiny inaccuracies.
	 * @param lowerThreshold25D is used to take also the heights of the
	 * involved points into account. For example, the points\n
	 * \n
	 * Point_2 a(0,0,0);\n
	 * Point_2 b(10,0,0);\n
	 * Point_2 c(10,10,0);\n
	 * Point_2 d(0,10,1000);\n
	 * \n
	 * can form the triangles (a,b,c) and (a,c,d) or the triangles (a,b,d)
	 * and (d,b,c) but (a,c) is obviousy the better diagonal because the
	 * points a,b,c share the same elevation while d is at z=1000.
	 * Technically spoken, the diagonal with the smaller 2.5D-length is
	 * applied if the both, the 2D error is below \p tolerance2D and the
	 * 2.5D error is above \p lowerThreshold25D. The 2.5D
	 * criterion has priority over the 2D criterion.
	 *
	 */
	void optimizeValleysAndRidges(double tolerance2D,double lowerThreshold25D);
#endif

	/**
	 * Unify Grid
	 *
	 * A Delaunay triangulation is not unique when when 2 or more triangles
	 * lie on a common circle. As a consequence the four corners of
	 * a rectangle can be triangulated in two different ways: Either the
	 * diagonal proceeds from the lower left to the upper right corner
	 * or it connects the other two corners. Both solutions are valid and
	 * an arbitrary one is applied when points are triangulated. To improve
	 * the repeatability and for reasons of visual appearance this method
	 * unifies such diagonals to point from the lower left to the upper
	 * right corner (or in horizontal direction).
	 *
	 * @param tolerance is 0 when only exact cases of more than 3 points
	 * on a common circumcircle shall be changed. But in practice input
	 * data can be disturbed by noise and tiny rounding errors such that
	 * grid points are not exactly on a grid. The numeric error is computed
	 * as \f$error=\frac{abs(diagonalA-diagonalB)}{max(diagonalA,diagonalB)}\f$.
	 * and \p tolerance is an upper threshold to allow modification despite
	 * such tiny inaccuracies. Use with caution, such flips break the
	 * empty circle property and this may or may not fit your setting.
	 */
	void unifyGrid(double tolerance);

	/// @private
	bool assignDt2(Dt2* pDt_);

	/** \brief Get the triangles of the zone.
	*
	* This command fetches the existing triangles of the zone.
	*
	* @note Fade_2D::void applyConstraintsAndZones() must be called after
	* the last insertion of points and constraints.
	*
	* @note that the lifetime of data from the Fade2D datastructures
	* does exceed the lifetime of the Fade2D object.
	*/
	void getTriangles(std::vector<Triangle2*>& vTriangles_) const;

	/** \brief Get the vertices of the zone.
	*
	*/
	void getVertices(std::vector<Point2*>& vVertices_) const;


	/** Statistics
	 *
	 * Prints statistics to stdout.
	 */
	void statistics(const char* s) const;


	/** \brief Get the associated constraint
	* @return a pointer to the ConstraintGraph2 object which defines the zone.@n
	* or NULL for ZL_RESULT-, ZL_GROW and ZL_GLOBAL_-zones.
	*/
	ConstraintGraph2* getConstraintGraph() const;

	/** \brief Get the number of triangles
	* @warning This method is fast but O(n), so don't call it frequently in a loop.
	*
	*/
	size_t getNumberOfTriangles() const;


	/** \brief Get the associated constraint graphs
	*
	*/
	void getConstraintGraphs(std::vector<ConstraintGraph2*>& vConstraintGraphs_) const;

	/// @private
	Dt2* getDelaunayTriangulation() const;

	/** \brief Get a the number of ConstraintGraph2 objects
	*
	* A Zone2 object might be defined by zero, one or more ConstraintGraph2 objects.
	*/
	size_t numberOfConstraintGraphs() const;

	/** \brief Development function
	 */
	void debug(const char* name="");

	/** \brief Compute the bounding box
	 */
	Bbox2 getBoundingBox() const;

	// Deprecated, replaced by getBorderEdges() but kept for
	// backwards-compatibility.
	/** @private
	 */
	void getBoundaryEdges(std::vector<Edge2>& vEdges) const;

	/** \brief Get connected components and their boundary polygons
	 *
	 * This method subdivides the zone into connected components. For
	 * each connected component it then returns a CompPolygon object
	 * consisting of the triangles, their outer boundary polygon and
	 * the hole polygons. Edges are represented by a triangle and an
	 * index and they are always counterclockwise (CCW) around the
	 * triangle belonging to the Zone2. Thus the outer boundary polygon
	 * is ccw-oriented while the polygons of inner holes are cw-oriented.
	 */
	void getComponentPolygons(std::vector<CompPolygon>& vCompPolygons) const;

	/** \brief Compute the boundary segments
	 *
	 * Outputs the boundary segments of the zone. These are ccw-oriented
	 * but not returned in any specific order.
	 */
	void getBoundarySegments(std::vector<Segment2>& vSegments) const;

	/** \brief Get 2D Area
	 *
	 * Returns the 2D area of the zone.
	 *
	 * \if SECTION_FADE25D
	 * Note: The getArea() method is deprecated and replaced by getArea2D()
	 * and getArea25D()
	 * \else
	 * Note: The getArea() method is deprecated and replaced by getArea2D()
	 * to keep the names consistent.
	 * \endif
	 */
	double getArea2D() const;


#if GEOM_PSEUDO3D==GEOM_TRUE
	/** \brief Get 2.5D Area
	 *
	 * Returns the 2.5D area of the zone.
	 *
	 * Note: The getArea() method is deprecated and replaced by getArea2D()
	 * and getArea25D()
	 */
	double getArea25D() const;
#endif

	/** \brief Get border edges
	 * @return: the CCW oriented border edges of the zone
	 */
	void getBorderEdges(std::vector<Edge2>& vBorderEdgesOut) const;

	/** \brief Write the zone to *.obj
	* Writes the triangles of the present Zone2 to an *.obj file (The
	* *.obj format represents a 3D scene).
	*
	* @param outFilename is the output filename
	*/
	void writeObj(const char* outFilename) const;

#ifndef SKIPTHREADS
#if GEOM_PSEUDO3D==GEOM_TRUE
/** \brief Write the zone to a *.ply file
*
* @param filename is the output filename
* @param bASCII specifies whether to write the *.ply in ASCII or binary format
 *
 * @note Method available for platforms with C++11
*/
	bool writePly(const char*  filename,bool bASCII=false) const;
/** \brief Write the zone to a *.ply file
*
* @param os is the output file
* @param bASCII specifies whether to write the *.ply in ASCII or binary format
 *
 * @note Method available for platforms with C++11
*/
	bool writePly(std::ostream& os,bool bASCII=false) const;
#endif
#endif



protected:
	Zone2& operator=(const Zone2&);
	// Optimization techniques
	/// @private
	void optMode_standard_sub(std::vector<Triangle2*>& vT,std::vector<Triangle2*>& vChangedT);
	/// @private
	void optMode_standard();
	/// @private
	double optMode_prioq(double noEdgeBelowDih,bool bWithProgress);
	/// @private
	void getEdgesForOptimization(double noEdgeBelowDegree,std::vector<Edge2>& vEdges);
	/// @private
	void optMode_simulatedAnnealing();
	/// @private
	void optMode_simulatedAnnealing_sub(std::vector<Edge2>& vUniqueEdges,double temperature);
	/// @private
	void removeConstraintEdges(std::vector<Edge2>& vEdges) const;
/// @private
	Zone2(const Zone2&);
/// @private
	void getTriangles_RESULT(std::vector<Triangle2*>& vTriangles) const;
/// @private
	void initWorkspace(bool bInside,std::set<std::pair<Point2*,Point2*> >& sNoGrowEdges,std::vector<Triangle2*>& vWorkspace) const;
/// @private
	void bfsFromWorkspace(std::vector<Triangle2*>& vWorkspace,std::set<std::pair<Point2*,Point2*> >& sNoGrowEdges,std::vector<Triangle2*>& vTriangles) const;
/// @private
	Zone2* ctbz_treatCC(std::vector<Triangle2*>& vOneCC);
	// Data
/// @private
	Dt2* pDt;
/// @private
	Progress* pZoneProgress;
/// @private
	ZoneShooter* pZoneShooter;
/// @private
	ZoneLocation zoneLoc;
	CLASS_DECLSPEC
	friend Zone2* zoneUnion(Zone2* pZone0,Zone2* pZone1);
	CLASS_DECLSPEC
	friend Zone2* zoneIntersection(Zone2* pZone0,Zone2* pZone1);
	CLASS_DECLSPEC
	friend Zone2* zoneDifference(Zone2* pZone0,Zone2* pZone1);
	CLASS_DECLSPEC
	friend Zone2* zoneSymmetricDifference(Zone2* pZone0,Zone2* pZone1);

	/** \brief Peel off border triangles (deprecated)
	*
	* This function is DEPRECATED but kept for backwards compatibility.
	* The new and better function is:
	* peelOffIf(Zone2* pZone, bool bAvoidSplit,PeelPredicateTS* pPredicate)
	*
	* @param pZone
	* @param pPredicate
	* @param bVerbose
	* @return a new zone containing a subset of the triangles of \p pZone or NULL if no triangles remain.
	*/
	CLASS_DECLSPEC
	friend Zone2* peelOffIf(Zone2* pZone, UserPredicateT* pPredicate,bool bVerbose); // Depricated!

	/** \brief Peel off border triangles
	 *
	 * @param pZone is the input zone
	 * @param bAvoidSplit if true, then the algorithm removes a
	 * triangle only if it does not break the zone into independent
	 * components.
	 * @param pPredicate is a user-defined predicate that decides
	 * if a triangle shall be removed.
	 *
	 * @return a new zone containing a subset of the triangles of \p pZone or NULL if no triangles remain.
	 */
	CLASS_DECLSPEC
	friend Zone2* peelOffIf(Zone2* pZone, bool bAvoidSplit,PeelPredicateTS* pPredicate);



private:
#ifndef __MINGW32__
#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable:4251)
#endif
#endif
	std::vector<ConstraintGraph2*> vConstraintGraphs;
	std::vector<Point2> vStartPoints;
	std::vector<Zone2*> vInputZones;

#ifndef __MINGW32__
#ifdef _WIN32
#pragma warning(pop)
#endif
#endif
};

// Free functions
/**
 * @brief Computes the union of two Zone2 objects.
 *
 * This function returns a pointer to a new Zone2 object representing the union
 * of the two input zones. The union is the area covered by either of the input
 * zones. The input zones must belong to the same Fade_2D object.
 *
 * @param pZone0 Pointer to the first Zone2 object.
 * @param pZone1 Pointer to the second Zone2 object.
 * @return Zone2* A pointer to the resulting Zone2 object representing the union.
 */
CLASS_DECLSPEC
Zone2* zoneUnion(Zone2* pZone0,Zone2* pZone1);


/**
 * @brief Computes the intersection of two Zone2 objects.
 *
 * This function returns a pointer to a new Zone2 object representing the intersection
 * of the two input zones. The intersection is the area covered by both of the input
 * zones. The input zones must belong to the same Fade_2D object.
 *
 * @param pZone0 Pointer to the first Zone2 object.
 * @param pZone1 Pointer to the second Zone2 object.
 * @return Zone2* A pointer to the resulting Zone2 object representing the intersection.
 * If the result is empty, the function returns `NULL`.
 */
CLASS_DECLSPEC
Zone2* zoneIntersection(Zone2* pZone0,Zone2* pZone1);

/**
 * @brief Computes the difference between two Zone2 objects.
 *
 * This function returns a pointer to a new Zone2 object representing the difference
 * between the two input zones. The difference is the area covered by the first zone
 * that is not covered by the second zone. The input zones must belong to the same
 * Fade_2D object.
 *
 * @param pZone0 Pointer to the first Zone2 object.
 * @param pZone1 Pointer to the second Zone2 object.
 * @return Zone2* A pointer to the resulting Zone2 object representing the difference.
 * If the result is empty, the function returns `NULL`.
 */
CLASS_DECLSPEC
Zone2* zoneDifference(Zone2* pZone0,Zone2* pZone1);

/**
 * @brief Computes the symmetric difference between two Zone2 objects.
 *
 * This function returns a pointer to a new Zone2 object representing the symmetric
 * difference between the two input zones. The symmetric difference is the area covered
 * by either of the input zones but not by both. The input zones must belong to the
 * same Fade_2D object.
 *
 * @param pZone0 Pointer to the first Zone2 object.
 * @param pZone1 Pointer to the second Zone2 object.
 * @return Zone2* A pointer to the resulting Zone2 object representing the symmetric difference.
 * If the result is empty, the function returns `NULL`.
 */
CLASS_DECLSPEC
Zone2* zoneSymmetricDifference(Zone2* pZone0,Zone2* pZone1);




} // (namespace)
