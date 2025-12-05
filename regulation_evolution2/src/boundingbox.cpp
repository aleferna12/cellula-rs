//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include <cmath>
#include "boundingbox.h"

BoundingBox::BoundingBox(vector<pair<double, double>> const &points)
    : BoundingBox((int) round(points[0].first),
                  (int) round(points[0].second),
                  (int) round(points[0].first),
                  (int) round(points[0].second)) {
    for (auto &p : points) {
        if (p.first < minx)
            minx = (int) floor(p.first);
        else if (p.first > maxx)
            maxx = (int) ceil(p.first);
        if (p.second < miny)
            miny = (int) floor(p.second);
        else if (p.second > maxy)
            maxy = (int) ceil(p.second);
    }
}

pair<int, int> BoundingBox::getOverlapLengths(const BoundingBox &box) const {
    int lenx = max(minx, box.minx) - min(maxx, box.maxx);
    int leny = max(miny, box.miny) - min(maxy, box.maxy);
    return {lenx, leny};
}

int BoundingBox::getOverlapArea(const BoundingBox &box) const {
    auto lens = getOverlapLengths(box);
    return max(0, lens.first * lens.second);
}

bool BoundingBox::overlaps(const BoundingBox &box) const {
    auto lens = getOverlapLengths(box);
    return lens.first > 0 and lens.second > 0;
}
