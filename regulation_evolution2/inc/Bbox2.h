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
/// \file Bbox2.h

#pragma once
#include "Segment2.h"
#include "common.h"
#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

class GeomTest; // FWD


/** \brief Bbox2 is an axis aligned 2D bounding box.
*
*/
class CLASS_DECLSPEC Bbox2
{
public:
/** \brief Constructor
 *
 * Minimum bounds are initialized to DBL_MAX.
 * Maximum bounds are initialized to -DBL_MAX.
 * Box is not valid yet
*/
	explicit Bbox2(GeomTest* pGeomTest_=NULL):
		minX(DBL_MAX),minY(DBL_MAX),
		maxX(-DBL_MAX),maxY(-DBL_MAX),
		maxAbsCoord(0),
		bValid(false),pGeomTest(pGeomTest_)
	{
	}
	~Bbox2();


/** \brief Reset the bounds
 *
 * Resets the bounding box to default values
*/
	void reset()
	{
		minX=DBL_MAX;
		minY=DBL_MAX;
		maxX=-DBL_MAX;
		maxY=-DBL_MAX;
		bValid=false;
	}


/** \brief Check if the bounds are valid
 *
 * The bounds are valid when at least one point has been added or
 * when set-methods have been used to set minX<=maxX and minY<=maxY
*/

	bool isValid() const
	{
		return minX<DBL_MAX;
	}

/** \brief Get corners
 *
 * Convenience function: Returns the 4 corners of the box
 */
	void getCorners(std::vector<Point2>& vBoxCorners) const;

/** \brief Get boundary
 *
 * Convenience function: Returns 4 border segments
 */
	void getBoundary(std::vector<Segment2>& vBoundary) const;

/** \brief Get offset corners
 *
 * Convenience function: Returns the 4 corners of an enlarged box. The
 * box es enlarged by \p offset in each direction
 */
	void getOffsetCorners(double offset,std::vector<Point2>& vBoxCorners) const;


/** \brief Check intersection
 *
 * Two valid bounding boxes intersect if they share at least one point in the XY plane.
 */
	bool doIntersect(const Bbox2& other) const;

/** \brief Add points
 *
 * Extends the 2D bounding box if required.
* \return true if the bounding box changes, false otherwise
*/
	bool add(std::vector<Point2*>::const_iterator start_it,std::vector<Point2*>::const_iterator end_it)
	{
		if(start_it==end_it) return false;
		if(bValid)
		{
			double oldMinX(minX),oldMinY(minY),oldMaxX(maxX),oldMaxY(maxY);
			for(;start_it!=end_it;++start_it) treatPointForValidBox(**start_it);
			updateMaxAbs();
			if(oldMinX!=minX || oldMinY!=minY || oldMaxX!=maxX || oldMaxY!=maxY ) return true;
			return false;
		}
		else
		{
			treatPointForInvalidBox(**start_it);
			++start_it;
			for(;start_it!=end_it;++start_it) treatPointForValidBox(**start_it);
			updateMaxAbs();
			return true;
		}

	}
/** \brief Add points
 *
 * Extends the 2D bounding box if required.
* \return true if the bounding box changes, false otherwise
*/

	bool add(std::vector<Point2>::const_iterator start_it,std::vector<Point2>::const_iterator end_it)
	{
		if(start_it==end_it) return false;
		if(bValid)
		{
			double oldMinX(minX),oldMinY(minY),oldMaxX(maxX),oldMaxY(maxY);
			for(;start_it!=end_it;++start_it) treatPointForValidBox(*start_it);
			updateMaxAbs();
			if(oldMinX!=minX || oldMinY!=minY || oldMaxX!=maxX || oldMaxY!=maxY ) return true;
			return false;
		}
		else
		{
			treatPointForInvalidBox(*start_it);
			++start_it;
			for(;start_it!=end_it;++start_it) treatPointForValidBox(*start_it);
			updateMaxAbs();
			return true;
		}
	}

/** \brief Add points
 *
 * Extends the 2D bounding box if required.
* \return true if the bounding box changes, false otherwise
*/

	bool add(size_t numPoints,double * coordinates)
	{
#if GEOM_PSEUDO3D==GEOM_TRUE
		const int NUMCOMPONENTS(3);
#else
		const int NUMCOMPONENTS(2);
#endif

		if(numPoints==0) return false;
		double oldMinX(minX),oldMinY(minY),oldMaxX(maxX),oldMaxY(maxY);
		double firstX(coordinates[0]);
		double firstY(coordinates[1]);
		if(firstX<minX) minX=firstX;
		if(firstX>maxX) maxX=firstX;
		if(firstY<minY) minY=firstY;
		if(firstY>maxY) maxY=firstY;

		for(size_t i=0;i<numPoints;++i)
		{
			double x(coordinates[NUMCOMPONENTS*i]);
			double y(coordinates[NUMCOMPONENTS*i+1]);
			if(x<minX) minX=x;
				else if(x>maxX) maxX=x;
			if(y<minY) minY=y;
				else if(y>maxY) maxY=y;
		}
		bValid=true;
		updateMaxAbs();
		if(oldMinX!=minX || oldMinY!=minY || oldMaxX!=maxX || oldMaxY!=maxY ) return true;
			else return false;
	}

/** \brief Add a point
 *
 * Extends the 2D bounding box if required.
* \return true if the bounding box changes, false otherwise
*/

	bool add(const Point2& p)
	{
		//std::cout<<"Add point: "<<p<<std::endl;
		if(bValid)
		{
			double oldMinX(minX),oldMinY(minY),oldMaxX(maxX),oldMaxY(maxY);
			treatPointForValidBox(p);
			updateMaxAbs();
			if(oldMinX!=minX || oldMinY!=minY || oldMaxX!=maxX || oldMaxY!=maxY ) return true;
				else return false;
		}
		else
		{
			treatPointForInvalidBox(p);
			updateMaxAbs();
			return true;
		}
	}
/** \brief Add another bounding box
 *
 * Extends the 2D bounding box if required.
* \return true if the bounding box changes, false otherwise
*/
	bool add(const Bbox2& other)
	{
		bool bRet(false);
		if(other.minX<minX) {minX=other.minX;bRet=true;}
		if(other.minY<minY) {minY=other.minY;bRet=true;}
		if(other.maxX>maxX) {maxX=other.maxX;bRet=true;}
		if(other.maxY>maxY) {maxY=other.maxY;bRet=true;}
		updateMaxAbs();
		return bRet;
	}

/** \brief Point-in-Box Test
 *
 * @return true if minX <= p.x() <=maxX and minY <= p.y() <=maxY or false otherwise.
*/
	bool isInBox(const Point2& p) const;

/** \brief Compute the 2D midpoint
 *
 */
	Point2 computeCenter() const;

/** \brief Add a bounding box
 *
 * Extends the 2D bounding box if required.
* \return the resulting bounding box
*/
	Bbox2 operator+(const Bbox2& b);


#if GEOM_PSEUDO3D==GEOM_TRUE
/** \brief Get the min point
 *
* \return the corner point with the minimum coordinates, the z-coordinate is set to 0
*/
#else
/** \brief Get the min point
 *
* \return the 2D corner point with the minimum coordinates
*/
#endif

	Point2 getMinPoint() const
	{
#if GEOM_PSEUDO3D==GEOM_TRUE
		return Point2(minX,minY,0);
#else
		return Point2(minX,minY);
#endif
	}


#if GEOM_PSEUDO3D==GEOM_TRUE
/** \brief Get the max point
 *
* \return the 2D corner point with the maximum coordinates, the z-coordinate is set to 0
*/
#else
/** \brief Get the max point
 *
* \return the 2D corner point with the maximum coordinates
*/
#endif

	Point2 getMaxPoint() const
	{
#if GEOM_PSEUDO3D==GEOM_TRUE
		return Point2(maxX,maxY,0);
#else
		return Point2(maxX,maxY);
#endif
	}

/** \brief Get minimum coordinate
*
* \return the smallest coordinate value, i.e. min(minX,minY)
*/

	double getMinCoord() const
	{
		return (std::min)(minX,minY);
	}
/** \brief Get maximum coordinate
 *
* \return the largest coordinate value, i.e. max(maxX,maxY)
*/

	double getMaxCoord() const
	{
		return (std::max)(maxX,maxY);
	}

/** \brief Get x-range
 *
* \return maxX-minX
*/

	double getRangeX() const
	{
		return maxX-minX;
	}
/** \brief Get y-range
 *
* \return maxY-minY
*/
	double getRangeY() const
	{
		return maxY-minY;
	}

/** \brief Get max range
*
* \return the largest range, i.e. max(getRangeX(),getRangeY())
*/
	double getMaxRange() const
	{
		double range0=getRangeX();
		double range1=getRangeY();
		if(range0>range1) return range0;
		return range1;
	}

/** \brief Get minX
 *
* \return minX
*/
	double get_minX() const {return minX;}
/** \brief Get minY
 *
* \return minY
*/
	double get_minY() const {return minY;}
/** \brief Get maxX
 *
* \return maxX
*/
	double get_maxX() const {return maxX;}
/** \brief Get maxY
 *
* \return maxY
*/
	double get_maxY() const {return maxY;}

/** \brief Get bounds
 *
*/
	void getBounds(double& minX_,double& maxX_,double& minY_,double& maxY_) const;


/** \brief Double the box
 *
* Changes the bounds such that the box grows in each direction by half
* the previous range
*/
	void doubleTheBox();

/** \brief Set minX
*/
	void setMinX(double val)
	{
		minX=val;
		updateMaxAbs();
		if(minX<=maxX && minY<=maxY) bValid=true;
	}
/** \brief Set maxX
*/
	void setMaxX(double val)
	{
		maxX=val;
		updateMaxAbs();
		if(minX<=maxX && minY<=maxY) bValid=true;
	}
/** \brief Set minY
*/
	void setMinY(double val)
	{
		minY=val;
		updateMaxAbs();
		if(minX<=maxX && minY<=maxY) bValid=true;
	}
/** \brief Set maxY
*/
	void setMaxY(double val)
	{
		maxY=val;
		updateMaxAbs();
		if(minX<=maxX && minY<=maxY) bValid=true;
	}

/** @brief Enlarge the x|y-ranges of a bounding box
 *
 * This function updates \p minX, \p minY, \p maxX and \p maxY
 * symmetrically to enlarge the \p x- and \p y-ranges of the bounding
 * box.
 *
 * @param factor is the factor by which the existing ranges \p rangeX
 * and \p rangeY shall grow.
 * @param bUseMaxRange specifies if a square bounding box
 * with side length \f$ factor * max(rangeX,rangeY) \f$ shall be made.
 * @param minRange is a lower bound on the new ranges and it can be
 * used to avoid a degenerate box when the pre-existing box was degenerate.
 */
	void enlargeRanges(double factor,bool bUseMaxRange,double minRange);

/** \brief Inflate if Degenerate
 *
 * When only one point has been added to Bbox2 or when all points
 * have the same x- and/or y- coordinates then Bbox2 is degenerate.
 * This is a valid state but sometimes undesireable. The present
 * method inflates the Bbox2 by adding /p val to maxX and/or maxY.
*/
	void inflateIfDegenerate(double val)
	{
		if(bValid)
		{
			if(minX==maxX) maxX+=val;
			if(minY==maxY) maxY+=val;
		}
		updateMaxAbs();
	}
protected:
	inline void treatPointForValidBox(const Point2& p)
	{
		double x,y;
		p.xy(x,y);
		if(x<minX) minX=x;
			else if(x>maxX) maxX=x;
		if(y<minY) minY=y;
			else if(y>maxY) maxY=y;
		updateMaxAbs();
	}
	inline void treatPointForInvalidBox(const Point2& p)
	{
		// Individual bounds may have been set already. Keep them!
		if(minX==DBL_MAX) minX=p.x();
		if(minY==DBL_MAX) minY=p.y();
		if(maxX==-DBL_MAX) maxX=p.x();
		if(maxY==-DBL_MAX) maxY=p.y();
		updateMaxAbs();
		bValid=true;
	}
	friend std::ostream &operator<<(std::ostream &stream, const Bbox2& pC);
protected:
	double minX,minY;
	double maxX,maxY;
	double maxAbsCoord;
	bool bValid;
	GeomTest* pGeomTest;
	bool updateMaxAbs();
	friend Bbox2 intersection(const Bbox2& a,const Bbox2& b);
};

/** @brief Compute the bounding box
 *
 * Computes the bounding box of points
 */
Bbox2 getBox(std::vector<Point2>& vP);
/** @brief Compute the bounding box
 *
 * Computes the bounding box of points
 */
Bbox2 getBox(std::vector<Point2*>& vP);
/** @brief Print the box
 *
 * Prints the box coordinates to \p stream
 */

inline Bbox2 intersection(const Bbox2& a, const Bbox2& b)
{
	Bbox2 ret(a);
	if(b.minX > ret.minX) ret.minX=b.minX;
	if(b.maxX < ret.maxX) ret.maxX=b.maxX;
	if(b.minY > ret.minY) ret.minY=b.minY;
	if(b.maxY < ret.maxY) ret.maxY=b.maxY;
	if( (ret.minX>ret.maxX) || (ret.minY>ret.maxY) )
	{
		ret.bValid=false;
		ret.minX=DBL_MAX;
		ret.minY=DBL_MAX;
		ret.maxX=-DBL_MAX;
		ret.maxY=-DBL_MAX;
	}
	else
	{
		ret.bValid = true;
	}
	return ret;
}


inline std::ostream &operator<<(std::ostream &stream, const Bbox2& pC)
{
	stream<<"Bbox2: xy("<<pC.minX<<","<<pC.minY<<") -> xy("<<pC.maxX<<","<<pC.maxY<<"), rangeX="<<pC.getRangeX()<<", rangeY="<<pC.getRangeY();
	return stream;
}



} // (namespace)
