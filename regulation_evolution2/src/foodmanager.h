//
// Created by aleferna on 31/10/23.
//

#ifndef REGULATION_EVOLUTION_FOODMANAGER_H
#define REGULATION_EVOLUTION_FOODMANAGER_H

#include <algorithm>
#include "voronoidiagram.h"
#include "foodpatch.h"
#include "intplane.h"

struct FoodPatchInfluence {
    BoundingBox bb;
    vector<int> neighbors;
};

class FoodManager {
public:
    FoodManager(int sizex,
                int sizey,
                int maxfoodpatches,
                int foodperspot,
                int foodpatcharea,
                int foodpatchperiod,
                int seasonduration,
                double seasonamplitude,
                double gradscale,
                double gradnoise);

    ~FoodManager();

    // Was too lazy to define copy constructors so now this is a singleton (which is prob. good anyway)
    FoodManager(FoodManager const &) = delete;

    FoodManager &operator =(FoodManager const &) = delete;

    int getLastAdded() const {
        return lastadded;
    }

    void setLastAdded(int value) {
        lastadded = value;
    }

    int getFoodLeft() const {
        int food = 0;
        for (auto &fp: foodpatches)
            food += fp.getFoodLeft();
        return food;
    }

    int getFoodPatchLength() const {
        return foodpatchlength;
    }

    int foodSigma(int x, int y) const {
        return foodplane->Sigma(x, y);
    }

    int chemSigma(int x, int y) const {
        return chemplane->Sigma(x, y);
    }

    const vector<int> &getActiveFoodPatches() const {
        return active_foodpatches;
    }

    FoodPatch const &getFoodPatch(int id) const {
        return foodpatches.at(id);
    }

    int numberFoodPatches() const {
        return (int) active_foodpatches.size();
    }

    tuple<int, int> minMaxChemSignal() const {
        return chemplane->getMinMax();
    }

    // Get distance of coordinate to the center of the closest food patch
    tuple<int, double> closestFoodPatch(int x, int y, vector<int> const &check) const;

    // ChemPlane in relation to the distance from a single peak F(x)
    double chemEquation(double dist_from_peak) const;

    void replenishFood(int time);

    int addFoodPatch(int centerx, int centery, int *sigmas = nullptr);

    int addRandomFoodPatch();

    void removeFoodPatch(int fp_id);

    double determineMinDist(int n) const;

    double distMostIsolatedPoint() const;

    void initSigmas(int fp_id, const int *sigmas = nullptr);

    int consumeFood(int fp_id, int gi, int gj);

    bool checkEmpty(int fp_id) const;

    void updateChemSignal(int fp_id, vector<int> const &check);

private:
    vector<FoodPatch> foodpatches{};
    VoronoiDiagram voronoidiagram{};
    // Current non-empty food patches, bookkeeping is done in UpdateVoronoi()
    vector<int> active_foodpatches{};
    IntPlane *foodplane;
    IntPlane *chemplane;

    int maxfoodpatches;
    int foodperspot;
    int foodpatchlength;
    int foodpatcharea;
    int foodpatchperiod;
    int seasonduration;
    double seasonamplitude;
    double gradscale;
    double gradnoise;
    vector<pair<int, int>> circle_pos;
    int lastadded = 0;

    // Determines how much food is in a specific position
    double chemAtPosition(int x, int y, vector<int> const &check) const;

    void updateActiveFoodPatches();

    bool isTimeForFood(int time) const;

    static vector<pair<int, int>> circlePositions(int diameter);
};


#endif //REGULATION_EVOLUTION_FOODMANAGER_H
