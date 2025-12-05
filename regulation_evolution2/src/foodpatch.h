//
// Created by aleferna on 31-10-2022.
//

#ifndef REGULATION_EVOLUTION_FOODPATCH_H
#define REGULATION_EVOLUTION_FOODPATCH_H

#include <vector>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include "boundingbox.h"

class FoodPatch {
public:
    //! Default constructor. Parameter sigmas can be optionally provided to point to an iterator
    //! that contains buffered foodSigma values.
    FoodPatch(int id, int x, int y, int length, int food_per_spot);

    FoodPatch(FoodPatch const &fp);

    FoodPatch &operator =(FoodPatch const &fp);

    ~FoodPatch();

    int getX() const {
        return centerx - length / 2;
    }

    int getY() const {
        return centery - length / 2;
    }

    int getCenterX() const {
        return centerx;
    }

    int getCenterY() const {
        return centery;
    }

    int getId() const {
        return id;
    }

    int getLength() const {
        return length;
    }

    int getFoodPerSpot() const {
        return food_per_spot;
    }

    int getFoodLeft() const {
        return food_left;
    }

    bool isEmpty() const {
        return empty;
    }

    void addFoodLeft(int food) {
        food_left += food;
    }

    void updateFoodLeft() {
        food_left = 0;
        for (int i = 0; i < length * length; ++i)
            food_left += sigma[i];
    }

    int getSigma(int i, int j) const {
        assertInBounds(i, j);
        return sigma[i * length + j];
    }

    vector<int> getSigmasAsVector() const {
        vector<int> v {};
        v.reserve(length * length);
        for (int i = 0; i < length * length; i ++)
                    v.push_back(sigma[i]);
        return v;
    }

    int getGlobalX(int i) const {
        return i + getX();
    }

    int getGlobalY(int j) const {
        return j + getY();
    }

    int getLocalX(int i) const {
        return i - getX();
    }

    int getLocalY(int j) const {
        return j - getY();
    }

private:
    friend class FoodManager;
    int id;
    int centerx;
    int centery;
    int length;
    BoundingBox bb;
    int food_per_spot;
    int food_left = 0;
    int *sigma;
    // Whether there is any food spots for this patch
    // If true we can remove FoodPatch from the ChemPlane
    // If false, removed is also false
    bool empty = false;
    // Whether FoodPatch has been removed from the ChemPlane
    // If true, we can recycle the position on fpatches
    // If true, empty is also true
    bool removed = false;

    void setSigma(int i, int j, int val) {
        assertInBounds(i, j);
        sigma[i * length + j] = val;
    }

    bool inBounds(int i, int j) const {
        return bb.inside(i, j);
    }

    void assertInBounds(int i, int j) const {
        if (!inBounds(i, j))
            throw out_of_range(
                    "Out of bounds access to FoodPatch sigma array at x" + to_string(i)
                    + " y" + to_string(j));
    }
};


#endif //REGULATION_EVOLUTION_FOODPATCH_H
