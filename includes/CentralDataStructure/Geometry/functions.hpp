#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_FUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_FUNCTIONS_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"

#include <array>
#include <vector>

namespace cds
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

    inline double squaredDistance(const Coordinate& a, const Coordinate& b)
    {
        return squaredLength(a - b);
    }

    inline double dotProduct(const Coordinate& a, const Coordinate& b)
    {
        auto dot = [&](int n)
        {
            return a.nth(n) * b.nth(n);
        };
        return dot(0) + dot(1) + dot(2);
    }

    inline bool withinDistance(double distance, const Coordinate& a, const Coordinate& b)
    {
        return squaredDistance(a, b) < distance * distance;
    }

    inline bool withinSphere(const Sphere& sphere, const Coordinate& point)
    {
        return withinDistance(sphere.radius, sphere.center, point);
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
    Coordinate coordinateMean(const std::vector<Coordinate*>& coords);
} // namespace cds
#endif
