//
// Created by aleferna on 07/12/23.
//

#include <Fade_2D.h>
#include "voronoidiagram.h"

using namespace GEOM_FADE2D;

void VoronoiDiagram::generate(const vector<pair<double, double>> &ps, const BoundingBox bb, const vector<int> &indexes) {
    dt.reset();
    cells.clear();

    // Scaling the box by a factor of 10 per Fade2 suggestion: https://www.geom.at/example12-voronoi-diagram/
    double bbx = bb.getMaxX() - bb.getMinX();
    double bby = bb.getMaxY() - bb.getMinY();
    vector<Point2> vPoints {
        {bbx / 2 * (1 - 10) - bb.getMinX(), bby / 2 * (1 - 10) - bb.getMinY()},
        {bbx / 2 * (1 - 10) - bb.getMinX(), bby / 2 * (1 + 10) - bb.getMinY()},
        {bbx / 2 * (1 + 10) - bb.getMinX(), bby / 2 * (1 - 10) - bb.getMinY()},
        {bbx / 2 * (1 + 10) - bb.getMinX(), bby / 2 * (1 + 10) - bb.getMinY()},
    };

    int c = 0;
    for (auto const &p : ps) {
        Point2 pp(p.first, p.second);
        pp.setCustomIndex(indexes[c++]);
        vPoints.push_back(pp);
    }
    dt.insert(vPoints);

    pVoro = dt.getVoronoiDiagram();
    if (!pVoro->isValid())
        throw runtime_error("Invalid voronoi diagram");

    vector<VoroCell2*> vCells;
    pVoro->getVoronoiCells(vCells);
    for (auto cell : vCells) {
        auto index = cell->getSite()->getCustomIndex();
        cell->setCustomCellIndex(index);
        cells.insert({index, cell});
    }
}
