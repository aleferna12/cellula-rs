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


/// @file PolygonTree.h
#pragma once

#include "common.h"
#include "Segment2.h"
#include "Edge2.h"
#include "VertexPair2.h"

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif

/**
 * @brief The PolygonTree class represents nested polygon layers as
 * a hierarchical structure
 *
 * This class organizes nested polygons in a tree-like format, where
 * each layer of the tree corresponds to a polygon layer, progressing
 * from the outermost to the innermost. The direct children of the root
 * node represent layer 0, which corresponds to the outermost polygons
 * (potentially multiple in case of disjoint components). The leaf nodes
 * of the tree correspond to the innermost polygons.
 */
class CLASS_DECLSPEC PolygonTree
{
public:
	/**
     * @brief Constructs a PolygonTree object
     *
     * @param layer_ The layer number this tree node represents.
     */
	explicit PolygonTree(int layer_);

    /**
     * @brief Destructor for the PolygonTree class.
     */
	~PolygonTree();

	/**
	 * @brief Retrieves the layer number of this PolygonTree node.
	 *
	 * @return The layer number associated with this node.
	 *
	 * The root node of the tree has a layer number of -1 and does not
	 * contain any segments. Its direct children, which have a layer number of 0,
	 * represent the outermost segments. There can be multiple layer-0 nodes
	 * for disjoint components of the polygon. The child nodes of the layer-0
	 * nodes correspond to layer 1, representing the next inner layer, and so on.
	 */
    int getLayer() const;

    /**
     * @brief Retrieves the child nodes of this PolygonTree node.
     *
     * @return A vector of pointers to the child PolygonTree nodes.
     */
	std::vector<PolygonTree*>& getChildren();

    /**
     * @brief Retrieves all child nodes recursively.
     *
     * This method populates the provided vector with all child nodes of this
     * PolygonTree node, traversing the tree recursively.
     *
     * @param vChildNodesRec A vector to be filled with the recursive child nodes.
     */
	void getChildrenRecursive(std::vector<PolygonTree*>& vChildNodesRec);
	/**
	* @brief Retrieves region-oriented segments.
	*
	* The boundary layer (layer 0) of a polygon defines a transition
	* from empty space to filled area (or from 'air' to 'material').
	* If a polygon contains a nested polygon, it defines a hole;
	* the hole polygons of layer 1 indicate a transition from material
	* back to air. More generally, considering multiple nested polygon layers
	* from outside to inside signifies alternating transitions between
	* air and material. Even layers represent borders from air to
	* material, while odd layers represent borders to holes ('air') within
	* that material.
	*
	* This method returns the segments of the current layer, ordered
	* and oriented counter-clockwise around the adjacent 'material'
	* region. Consequently, segments from even polygon layers are
	* oriented counter-clockwise, while those from odd layers are
	* oriented clockwise.
	*
	* @param vSegments A vector to be filled with the region-oriented segments.
	*/
	void getSegments_regionOriented(std::vector<Segment2>& vSegments) const;

    /**
     * @brief Retrieves counter-clockwise ordered and oriented segments.
     *
     * This method populates the provided vector with segments ordered
     * and oriented in a counter-clockwise direction.
     *
     * @param vSegments A vector to be filled with CCW-oriented segments.
     */
	void getSegments_CCW(std::vector<Segment2>& vSegments) const;

    /**
     * @brief Retrieves clockwise ordered and oriented segments.
     *
     * This method fills the provided vector with segments ordered
     * and oriented in a clockwise direction.
     *
     * @param vSegments A vector to be filled with CW-oriented segments.
     */
	void getSegments_CW(std::vector<Segment2>& vSegments) const;

	/// @private
	void setGeometrySortedPairs(std::vector<VertexPair2>& vGeometrySortedCW);
	/// @private
	void getPointerSortedPairs(std::vector<VertexPair2>& vFromCW_sorted) const;
private:
	/// @private
	//	const std::vector<VertexPair2>& getGeometrySortedPairsCW();
	explicit PolygonTree(const PolygonTree&);
	PolygonTree& operator=(const PolygonTree&);
	// DATA
	struct Impl;
	Impl* pImpl;///< Pointer to the implementation details of the PolygonTree.

};



} // (namespace)
