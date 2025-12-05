//
// Created by aleferna on 07/12/23.
//

#ifndef REGULATION_EVOLUTION_VORONOIDIAGRAM_H
#define REGULATION_EVOLUTION_VORONOIDIAGRAM_H

#include <cfloat>
#include <utility>
#include <vector>
#include <unordered_map>
#include <Fade_2D.h>
#include "boundingbox.h"

using namespace GEOM_FADE2D;

class VoronoiDiagram {
public:
    void generate(const vector<pair<double, double>> &ps, BoundingBox bb, const vector<int> &indexes);

    int size() {
        return (int) cells.size();
    }

    vector<pair<double, double>> getCellBoundaries(int id) {
        vector<Point2> vPoints;
        cells.at(id)->getBoundaryPoints(vPoints);
        vector<pair<double, double>> ret;
        for (auto const &p : vPoints) {
            ret.emplace_back(p.x(), p.y());
        }
        return ret;
    }

    int locateCell(double x, double y) {
        throwNull();
        auto cell = pVoro->locateVoronoiCell(Point2(x, y));
        return cell->getCustomCellIndex();
    }

    //! Gets all neighbours that are not boundary cells or itself
    vector<int> findNeighbours(int id) {
        vector<Point2 *> vNeighs;
        cells.at(id)->getNeighborSites(vNeighs);

        vector<int> neighs;
        neighs.reserve(vNeighs.size());
        for (auto const p : vNeighs) {
            auto index = p->getCustomIndex();
            if (index > -1 && index != id)
                neighs.push_back(index);
        }
        return neighs;
    }

private:
    void throwNull() {
        if (pVoro == nullptr)
            throw runtime_error("Uninitialized diagram, use 'generate' first");
    }

    Fade_2D dt{1};
    Voronoi2 *pVoro;
    unordered_map<int, VoroCell2 const *> cells = {};
    unordered_map<int, vector<pair<double, double>>> vertices = {};
    unordered_map<int, vector<int>> neighbours = {};
};


#endif //REGULATION_EVOLUTION_VORONOIDIAGRAM_H
