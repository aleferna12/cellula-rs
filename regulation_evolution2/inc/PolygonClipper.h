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


#pragma once
#include <vector>
#include "Segment2.h"
#include "PolygonTree.h"
#include "Visualizer2.h"
#include "common.h"

/// \file PolygonClipper.h

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif
class PolygonClipperImpl; // FWD
class Zone2; // FWD


/**
 * @brief Limits an input polygon to a specified zone.
 *
 * This function trims an input polygon such that its area is constrained to the given zone.
 * Optionally, the input segments can be draped onto the zone.
 *
 * @param pBaseZoneIn The zone that will be used to clip the polygon.
 * @param vPolygonInput A vector of unordered and unoriented Segment2s. The segments
 * may self-intersect, and the polygon can have holes. All such issues are automatically repaired.
 * @param bWithDrape Specifies whether the output polygon should be draped onto the zone.
 * @param vPolygonSegmentsOut A vector to store the resulting polygon segments. The output
 * segments are free of self-intersections and are oriented counterclockwise around the polygon area.
 */
bool CLASS_DECLSPEC clipPolygon(Zone2* pBaseZoneIn,std::vector<Segment2>& vPolygonInput,bool bWithDrape,std::vector<Segment2>& vPolygonSegmentsOut);

/**
 * @brief The PolygonClipper class handles polygon repair operations
 *
 * This class takes a set of polygon segments. After optional removement
 * of ultra-short segments, it resolves self-intersections and provides
 * methods to retrieve the layers of the repaired polygon.
 */
class CLASS_DECLSPEC PolygonClipper
{
public:

    /**
     * @brief Constructs a PolygonClipper object.
     *
     * Initializes the clipper with a set of segments representing polygon edges
     * and a collapse distance used to eliminate near-degenerate edges.
     *
     * @param vSegments A vector containing possibly intersecting polygon edges with no specific order or orientation.
     * @param collapseDist The distance threshold to collapse short edges. If this value is 0, short
     * edges are not collapsed. If this value is negative, a distance above the numeric uncertainty
     * (determined from the coordinates) is used.
     */
	PolygonClipper(const std::vector<Segment2>& vSegments, double collapseDist);
    /**
     * @brief Destructor for the PolygonClipper class.
     */
	~PolygonClipper();

	/**
	* @brief Visualizes the polygon areas as a .PDF or PostScript file.
	*
	* @param pVis A pointer to a Visualizer2 object.
	* @param matColor The color used for the 'material' regions.
	* @param airColor The color used for the 'air' regions (default is white and fully transparent).
	*/
	 void show(Visualizer2* pVis,const Color& matColor=Color(CYELLOW,0,true),const Color airColor=Color(CWHITE,0,false)) const;

	/**
	* @brief Visualizes the polygon areas as a .PDF or PostScript file.
	*
	* @param name The output filename.
	* @param matColor The color used for the 'material' regions.
	* @param airColor The color used for the 'air' regions (default is white and fully transparent).
	*/
	void show(const std::string& name,const Color& matColor=Color(CYELLOW,0,true),const Color airColor=Color(CWHITE,0,false)) const;


	/**
	* @brief Returnes the PolygonTree structure.
	*
	* This method returns a pointer to a `PolygonTree` object. This object
	* represents the root node of the hierarchical PolygonTree structure for
	* the repaired polygon. The structure can be analyzed layer by layer.
	*
	* @return A pointer to the root node of the PolygonTree.
	*/
	PolygonTree* getPolygonTree();
	/**
	* @brief Retrieves the outermost segments of the polygon in counter-clockwise direction.
	*
	* This method fills the provided vector with the polygon's outer boundary
	* segments, ordered and oriented in counter-clockwise (CCW) direction around
	* the enclosed area. Segments inside this outer polygon are ignored.
	*
	* @note The returned polygon is "traversed from the outside", meaning the largest
	* possible polygon is returned, rather than splitting a non-simple polygon
	* into multiple parts. However, the input polygon may consist of multiple
	* connected components, and in this case, more than one polygon is stored in
	* the output vector. Due to the counterclockwise orientation of the segments,
	* this output is still suitable as input for Fade_2D::createZone(ConstraintGraph2*,ZoneLocation,bool).
	*
	* @note If you need individual polygons for each connected component, you
	* may retrieve the root node of the PolygonTree using getPolygonTree(), and
	* then query its first-layer children (layer 0), which represent the individual
	* connected components of the outermost layer.
	*
	* @param vOuterSegments_CCW A vector to be filled with the CCW-oriented outer boundary segments.
	*/
	void getOuterSegments_CCW(std::vector<Segment2>& vOuterSegments_CCW) const;
	/**
	* @brief Retrieves the inner and outer polygon segments oriented by region.
	*
	* This method returns all layers of the polygon, oriented in a
	* counter-clockwise direction around their adjacent inside-regions.
	*
	* @param vSegments_regionOriented A vector to be filled with the region-oriented segments.
	*
	* **Even-odd rule:** The outermost layer (layer 0) defines a transition from empty
	* space to filled area (from 'air' to 'material'), while the subsequent
	* layer (layer 1) indicates holes (a transition back to 'air').
	* More generally, as we move from outside to inside, the layers
	* signify alternating transitions between air and material. Even
	* layers represent borders from air to material, while odd layers
	* represent borders to holes within that material. Consequently,
	* even layers are oriented counter-clockwise, while odd layers are
	* clockwise.
	*
	* The polygons are traversed to form the largest possible shapes
	* instead of splitting non-simple polygons into multiple parts.
	* All polygons are stored in the same output vector. Due to the
	* counter-clockwise orientation around the 'filled' sections, this
	* output is suitable as input for Fade_2D::createZone(ConstraintGraph2*,
	* ZoneLocation, bool).
	*
	* @note If you need individual polygons, retrieve the root node
	* of the PolygonTree using getPolygonTree() and query it
	* layer-by-layer. Each layer corresponds to one layer of the
	* polygon, from outside to inside.
	*/
	void getSegments_regionOriented(std::vector<Segment2>& vSegments_regionOriented) const;

private:
	/// @private
	PolygonClipper(const PolygonClipper &);
	/// @private
	PolygonClipper &operator=(const PolygonClipper &);


	// Data
	PolygonClipperImpl* pImpl;
};





} // Namespace

