#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMCOORDINATES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMCOORDINATES_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/references.hpp"

#include <vector>

namespace cds
{
    Sphere coordinateWithRadius(Atom* atom);
    std::vector<Coordinate> atomCoordinates(const std::vector<Atom*>& atoms);
    std::vector<CoordinateReference> atomCoordinateReferences(std::vector<Atom*>& atoms);
    std::vector<Sphere> atomCoordinatesWithRadii(const std::vector<Atom*>& atoms);
} // namespace cds
#endif