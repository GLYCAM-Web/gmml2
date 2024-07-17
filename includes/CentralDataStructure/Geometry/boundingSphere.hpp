#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"

#include <array>
#include <vector>

namespace cds
{
    struct Sphere
    {
        double radius;
        Coordinate center;
    };

    inline bool withinSphere(const Sphere& sphere, const Coordinate& point)
    {
        return withinDistance(sphere.radius, sphere.center, point);
    }

    inline bool spheresOverlap(double tolerance, const Sphere& a, const Sphere& b)
    {
        return withinDistance(a.radius + b.radius - tolerance, a.center, b.center);
    }

    Sphere boundingSphere(const std::vector<Coordinate>& points);
    Sphere boundingSphere(const std::vector<Sphere>& spheres);
} // namespace cds

#endif
