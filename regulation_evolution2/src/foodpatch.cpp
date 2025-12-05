//
// Created by aleferna on 31-10-2022.
//

#include <algorithm>
#include "foodpatch.h"
#include "misc.h"

FoodPatch::FoodPatch(
        int id,
        int x,
        int y,
        int length,
        int food_per_spot
        ) : id(id),
            centerx(x),
            centery(y),
            length(length),
            bb(0, 0, length - 1, length - 1),
            food_per_spot(food_per_spot) {
    sigma = new int[length * length]{};
}

FoodPatch::FoodPatch(
        const FoodPatch &fp
        ) : id(fp.id),
            centerx(fp.centerx),
            centery(fp.centery),
            length(fp.length),
            bb(fp.bb),
            food_per_spot(fp.food_per_spot),
            food_left(fp.food_left),
            empty(fp.empty),
            removed(fp.removed) {
    sigma = new int[length * length];
    copy(fp.sigma, fp.sigma + fp.length * fp.length, sigma);
}

FoodPatch &FoodPatch::operator=(const FoodPatch &fp) {
    if (this == &fp)
        return *this;
    if (length != fp.length)
        throw runtime_error("tried assigning FoodPatches of different length");

    id = fp.id;
    centerx = fp.centerx;
    centery = fp.centery;
    length = fp.length;
    food_per_spot = fp.food_per_spot;
    food_left = fp.food_left;
    empty = fp.empty;
    removed = fp.removed;

    copy(fp.sigma, fp.sigma + fp.length * fp.length, sigma);

    return *this;
}

FoodPatch::~FoodPatch() {
    delete[] sigma;
}
