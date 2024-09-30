#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_BOUNDINGSPHERE_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"

#include <vector>

namespace cds
{
    Sphere boundingSphereIncluding(Sphere sphere, const Sphere include);
    Sphere boundingSphere(const std::vector<Sphere>& spheres);
} // namespace cds

#endif
