#ifndef INCLUDES_GEOMETRY_BOUNDINGSPHERE_HPP
#define INCLUDES_GEOMETRY_BOUNDINGSPHERE_HPP

#include "include/geometry/geometryTypes.hpp"

#include <vector>

namespace gmml
{
    Sphere boundingSphereIncluding(Sphere sphere, const Sphere& include);
    Sphere boundingSphereCenteredOnLine(const Sphere& sphere, const Coordinate& point1, const Coordinate& point2);
    Sphere boundingSphere(const std::vector<Sphere>& spheres);
} // namespace gmml

#endif
