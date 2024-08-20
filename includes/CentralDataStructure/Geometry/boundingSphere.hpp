#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"

#include <vector>

namespace cds
{
    Sphere boundingSphere(const std::vector<Coordinate>& points);
    Sphere boundingSphere(const std::vector<Sphere>& spheres);
} // namespace cds

#endif
