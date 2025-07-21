#include "include/geometry/geometryFunctions.hpp"

#include "include/geometry/geometryTypes.hpp"

#include <cmath>
#include <vector>

namespace gmml
{
    double length(const Coordinate& a) { return std::sqrt(squaredLength(a)); }

    double distance(const Coordinate& a, const Coordinate& b) { return std::sqrt(squaredDistance(a, b)); }

    Coordinate scaleBy(double factor, const Coordinate& a)
    {
        auto scaled = [&](int n) { return factor * a.nth(n); };
        return {scaled(0), scaled(1), scaled(2)};
    }

    Coordinate normal(const Coordinate& a) { return scaleBy(1.0 / length(a), a); }

    Coordinate crossProduct(const Coordinate& a, const Coordinate& b)
    {
        auto cross = [&](int n, int k) { return a.nth(n) * b.nth(k) - a.nth(k) * b.nth(n); };
        return {cross(1, 2), cross(2, 0), cross(0, 1)};
    }

    Coordinate projection(const Coordinate& a, const Coordinate& b)
    {
        return scaleBy(dotProduct(a, b) / dotProduct(b, b), b);
    }

    Coordinate coordinateMean(const std::vector<Coordinate>& coords)
    {
        Coordinate center(0.0, 0.0, 0.0);
        for (auto& coord : coords)
        {
            center = center + coord;
        }
        return scaleBy(1.0 / coords.size(), center);
    }
} // namespace gmml
