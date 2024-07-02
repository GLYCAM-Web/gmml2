#ifndef INCLUDES_CENTRALDATASTRUCTURE_BOUNDINGSPHERE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_BOUNDINGSPHERE_HPP

#include "includes/CentralDataStructure/coordinate.hpp"

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

    Sphere boundingSphere(const std::vector<Coordinate>& points);
} // namespace cds

#endif
