#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <vector>

namespace cds
{
    Sphere boundingSphereIncluding(Sphere sphere, const Sphere& include);
    Sphere boundingSphereCenteredOnLine(const Sphere& sphere, const Coordinate& point1, const Coordinate& point2);
    Sphere boundingSphere(const std::vector<Sphere>& spheres);
} // namespace cds

#endif
