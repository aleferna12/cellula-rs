//
// Created by aleferna on 31/10/23.
//

#include <fstream>
#include "foodmanager.h"
#include "random.h"
#include "misc.h"

FoodManager::FoodManager(
        int sizex,
        int sizey,
        int maxfoodpatches,
        int foodperspot,
        int foodpatcharea,
        int foodpatchperiod,
        int seasonduration,
        double seasonamplitude,
        double gradscale,
        double gradnoise
        ) : maxfoodpatches(maxfoodpatches),
            foodperspot(foodperspot),
            foodpatcharea(foodpatcharea),
            foodpatchperiod(foodpatchperiod),
            seasonduration(seasonduration),
            seasonamplitude(seasonamplitude),
            gradscale(gradscale),
            gradnoise(gradnoise) {
    foodpatchlength = (int) ceil(sqrt(foodpatcharea / M_PI) * 2 + 2);
    foodpatchlength += foodpatchlength % 2 == 0;
    circle_pos = circlePositions(1024);
    foodplane = new IntPlane(sizex, sizey, -1);
    chemplane = new IntPlane(sizex, sizey, 0);
}

FoodManager::~FoodManager() {
    delete foodplane;
    delete chemplane;
}

tuple<int, double> FoodManager::closestFoodPatch(int x, int y, vector<int> const &check) const {
    double mindist_sq = foodplane->SizeX() * foodplane->SizeX() + foodplane->SizeY() * foodplane->SizeY();
    int res_id = -1;

    for (auto fp_id : check) {
        const auto &fp = foodpatches[fp_id];
        if (not fp.empty) {
            double dx = fp.getCenterX() - x;
            double dy = fp.getCenterY() - y;
            double dist_sq = dx * dx + dy * dy;
            if (dist_sq < mindist_sq) {
                mindist_sq = dist_sq;
                res_id = fp.getId();
            }
        }
    }
    return make_tuple(res_id, sqrt(mindist_sq));
}

// Alternatively, we could use the most isolated point to know where to put next peak at each iteration
// I think that doing this would be worse, as it is more computationally intensive and probably will tend to accumulate
// fpatches in the corners (?), which may be problematic for small grad_sources numbers
double FoodManager::determineMinDist(int n) const {
    double ratio = (foodplane->SizeY() - 2) / (double) (foodplane->SizeX() - 2);
    // ratio * sepx = sepy
    // sepx = sepy / ratio
    // (sepx + 1) * (sepy + 1) = grad_sources - 1
    // ratio * pow(sepx, 2) + sepx * (1 + ratio) + 2 - grad_sources = 0
    // Do the same for sepy and solve quadradic equations
    double sepx = solveQuadradic(ratio, 1 + ratio, 2 - n);
    double sepy = solveQuadradic(1 / ratio, 1 + 1 / ratio, 2 - n);
    double mindistx = (foodplane->SizeX() - 2) / (sepx * 2 + 2);
    double mindisty = (foodplane->SizeY() - 2) / (sepy * 2 + 2);
    return sqrt(mindistx * mindistx + mindisty * mindisty);
}

// TODO: This can be optimized to only search the voronoi cell it has to
double FoodManager::distMostIsolatedPoint() const {
    double dist = 0;
    for (int i = 1; i < foodplane->SizeX() - 1; i++)
        for (int j = 1; j < foodplane->SizeY() - 1; j++) {
            double closest_dist = get<1>(closestFoodPatch(i, j, active_foodpatches));
            if (closest_dist > dist) {
                dist = closest_dist;
            }
        }
    return dist;
}

double FoodManager::chemEquation(double dist_from_peak) const {
    // TODO: This is very slow because it computes sqrt twice for every position...
    // 1 is there to protect calculations from division by 0 errors
    return 1 + gradscale * foodplane->getDiagonal() / 100 * (1 - dist_from_peak / foodplane->getDiagonal());
}

double FoodManager::chemAtPosition(int x, int y, vector<int> const &check) const {
    double dist_from_peak = get<1>(closestFoodPatch(x, y, check));
    return chemEquation(dist_from_peak);
}

void FoodManager::removeFoodPatch(int fp_id) {
    auto check = voronoidiagram.findNeighbours(fp_id);
    updateChemSignal(fp_id, check);
    updateActiveFoodPatches();
    foodpatches[fp_id].removed = true;
}

int FoodManager::addFoodPatch(int centerx, int centery, int *sigmas) {
    // Center needs to be inside lattice because otherwise voronoi generation breaks
    if (centerx < 0 || centerx > foodplane->SizeX() - 1 || centery < 0 || centery > foodplane->SizeY() - 1)
        throw runtime_error("Tried to create food patch out of lattice limits");

    int fp_id = (int) foodpatches.size();
    for (auto &fp: foodpatches) {
        if (fp.removed) {
            fp_id = fp.getId();
            break;
        }
    }

    auto fp = FoodPatch(fp_id, centerx, centery, foodpatchlength, foodperspot);
    if (fp_id == (int) foodpatches.size())
        foodpatches.push_back(fp);
    else
        foodpatches.at(fp_id) = fp;
    initSigmas(fp_id, sigmas);

    if (!foodpatches.at(fp_id).empty) {
        updateActiveFoodPatches();
        updateChemSignal(fp_id, vector<int> {fp_id});
    }

    return fp_id;
}

// Adds a new FoodPatch at a semi-random position (still takes into account mindist)
int FoodManager::addRandomFoodPatch() {
    int centerx = 0, centery = 0;
    if (foodpatches.empty()) {
        return addFoodPatch(
            5 + (int) RandomNumber(foodplane->SizeX() - 10),
            5 + (int) RandomNumber(foodplane->SizeY() - 10)
        );
    }
    double dist = 0;
    double mindist = determineMinDist(int(foodpatches.size()) + 1);
    while (dist < mindist) {
        centerx = 5 + (int) RandomNumber(foodplane->SizeX() - 10);
        centery = 5 + (int) RandomNumber(foodplane->SizeY() - 10);
        dist = get<1>(closestFoodPatch(centerx, centery, active_foodpatches));
    }
    return addFoodPatch(centerx, centery);
}

// Initialize the FoodPatch on both the FoodPlane of dish and its internal plane
void FoodManager::initSigmas(int fp_id, const int *sigmas) {
    auto &fp = foodpatches.at(fp_id);
    int cx = fp.getCenterX();
    int cy = fp.getCenterY();
    int area = foodpatcharea;
    BoundingBox bb(2, 2, foodplane->SizeX() - 3, foodplane->SizeY() - 3);
    for (const auto &p : circle_pos) {
        int gi = p.first + cx;
        int gj = p.second + cy;

        if (!bb.inside(gi, gj))
            continue;

        int i = fp.getLocalX(gi);
        int j = fp.getLocalY(gj);

        if (!fp.inBounds(i, j)) {
            break;  // All positions that could be initialized were
        }

        int value;
        if (sigmas == nullptr)
            value = fp.getFoodPerSpot();
        else
            value = sigmas[i * fp.getLength() + j];
        if ((area > 0 || sigmas != nullptr) && value > 0 && foodplane->Sigma(gi, gj) == -1) {
            foodplane->setSigma(gi, gj, fp_id);
            fp.setSigma(i, j, value);
            fp.addFoodLeft(value);
            area--;
        }
    }
    // Just to be sure
    if (fp.getFoodLeft() == 0) {
        fp.empty = true;
        cerr << "Warning: tried to initialize FoodPatch at an invalid position: "
        << fp.getX() << ", " << fp.getY() << "\n";
    }
}

// Valid check to guarantee that there has been no miscount of food_left
bool FoodManager::checkEmpty(int fp_id) const {
    auto &fp = foodpatches.at(fp_id);
    int minx = max(2, fp.getX());
    int miny = max(2, fp.getY());
    int maxx = min(foodplane->SizeX() - 2, fp.getX() + fp.getLength());
    int maxy = min(foodplane->SizeY() - 2, fp.getY() + fp.getLength());
    for (int gi = minx; gi < maxx; ++gi) {
        for (int gj = miny; gj < maxy; ++gj) {
            if (foodplane->Sigma(gi, gj) == fp_id) {
                throw runtime_error("FoodPatch has no food left but is not empty");
            }
        }
    }
    return true;
}

// TODO: this should throw if food cant be consumed
int FoodManager::consumeFood(int fp_id, int gi, int gj) {
    if (foodplane->Sigma(gi, gj) != fp_id)
        throw runtime_error("consumeFood called with incorrect FoodPatch id");

    auto &fp = foodpatches.at(fp_id);

    int i = fp.getLocalX(gi);
    int j = fp.getLocalY(gj);
    int sigma_at = fp.getSigma(i, j);
    if (sigma_at > 0) {
        fp.setSigma(i, j, sigma_at - 1);
        fp.addFoodLeft(-1);
        if (sigma_at == 1) {
            foodplane->setSigma(gi, gj, -1);
            if (fp.getFoodLeft() == 0 and not fp.empty) {
                checkEmpty(fp_id);
                fp.empty = true;
                removeFoodPatch(fp.getId());
            }

        }
        return 1;
    }
    return 0;
}

bool FoodManager::isTimeForFood(int time) const {
    double fpperiod = sin(time * 2 * M_PI / seasonduration) * seasonamplitude + foodpatchperiod;
    return time - lastadded > fpperiod;
}

void FoodManager::replenishFood(int time) {
    // Try to add fpatch
    if (numberFoodPatches() < maxfoodpatches and isTimeForFood(time)) {
        addRandomFoodPatch();
        lastadded = time;
    }
}

void FoodManager::updateActiveFoodPatches() {
    active_foodpatches.clear();
    vector<pair<double, double>> points {};
    for (auto &fp : foodpatches) {
        if (fp.empty)
            continue;
        active_foodpatches.push_back(fp.getId());
        points.emplace_back(fp.getCenterX(), fp.getCenterY());
    }
    voronoidiagram.generate(
        points,
        BoundingBox(
            0,
            0,
            foodplane->SizeX() - 1,
            foodplane->SizeY() - 1
        ),
        active_foodpatches
    );
}

void FoodManager::updateChemSignal(int fp_id, vector<int> const &check) {
    auto bounds = voronoidiagram.getCellBoundaries(fp_id);
    auto area = BoundingBox(bounds);
    area.setMinX(max(0, area.getMinX()));
    area.setMinY(max(0, area.getMinY()));
    area.setMaxX(min(foodplane->SizeX() - 1, area.getMaxX()));
    area.setMaxY(min(foodplane->SizeY() - 1, area.getMaxY()));
    for (int i = area.getMinX(); i <= area.getMaxX(); i++)
        for (int j = area.getMinY(); j <= area.getMaxY(); j++) {
            if (voronoidiagram.locateCell(i, j) != fp_id)
                continue;

            double dfood = chemAtPosition(i, j, check);
            int local_maxfood = (int) dfood;
            chemplane->setSigma(i, j, local_maxfood);
            if (RANDOM() < dfood - local_maxfood) local_maxfood++;
            if (RANDOM() < gradnoise)
                chemplane->setSigma(i, j, local_maxfood);
        }
}

vector<pair<int, int>> FoodManager::circlePositions(int diameter) {
    int asize = diameter + (diameter % 2 == 0);
    double center = asize / 2.0;
    vector<pair<int, int>> pos;
    vector<double> dists;
    for (int i = 0; i < asize; i++) {
        for (int j = 0; j < asize; j++) {
            double dis = dist(i + 0.5, j + 0.5, center, center);
            dists.push_back(dis);
            if (dis > diameter / 2.0)
                continue;
            pos.emplace_back(i, j);
        }
    }
    sort(pos.begin(), pos.end(), [&dists, &asize](const pair<int, int> &p1, const pair<int, int> &p2) {
        double d1 = dists.at(p1.first * asize + p1.second);
        double d2 = dists.at(p2.first * asize + p2.second);
        if (d1 == d2) {
            return atan2(p1.first, p1.second) < atan2(p2.first, p2.second);
        }
        return d1 < d2;
    });
    for (auto &p : pos) {
        p.first -= diameter / 2;
        p.second -= diameter / 2;
    }
    return pos;
}
