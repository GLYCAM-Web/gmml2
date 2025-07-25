#ifndef INCLUDE_GEOMETRY_GEOMETRYFUNCTIONS_HPP
#define INCLUDE_GEOMETRY_GEOMETRYFUNCTIONS_HPP

#include "include/geometry/geometryTypes.hpp"

#include <array>
#include <vector>

namespace gmml
{
    inline double squaredLength(const Coordinate& a)
    {
        auto sq = [&](int n)
        {
            double d = a.nth(n);
            return d * d;
        };
        return sq(0) + sq(1) + sq(2);
    }

    inline double squaredDistance(const Coordinate& a, const Coordinate& b) { return squaredLength(a - b); }

    inline double dotProduct(const Coordinate& a, const Coordinate& b)
    {
        auto dot = [&](int n) { return a.nth(n) * b.nth(n); };
        return dot(0) + dot(1) + dot(2);
    }

    inline bool withinDistance(double distance, const Coordinate& a, const Coordinate& b)
    {
        double d = std::max(0.0, distance);
        return squaredDistance(a, b) < d * d;
    }

    inline bool spheresOverlap(double tolerance, const Sphere& a, const Sphere& b)
    {
        return withinDistance(a.radius + b.radius - tolerance, a.center, b.center);
    }

    double length(const Coordinate& a);
    double distance(const Coordinate& a, const Coordinate& b);
    Coordinate scaleBy(double factor, const Coordinate& a);
    Coordinate normal(const Coordinate& a);
    Coordinate crossProduct(const Coordinate& a, const Coordinate& b);
    // finds the point closest to 'a' on the line defined by 'b' and (0,0,0)
    Coordinate projection(const Coordinate& a, const Coordinate& b);
    Coordinate coordinateMean(const std::vector<Coordinate>& coords);
} // namespace gmml

#endif
