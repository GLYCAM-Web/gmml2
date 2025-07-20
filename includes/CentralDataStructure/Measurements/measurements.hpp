#ifndef INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <array>
#include <vector>

namespace cds
{
    Coordinate coordinateOppositeToNeighborAverage(
        const Coordinate& centralCoord, const std::vector<Coordinate>& neighbors, const double distance);

    Coordinate calculateCoordinateFromInternalCoords(
        const Coordinate& a,
        const Coordinate& b,
        const Coordinate& c,
        double angle_Degrees,
        double dihedral_Degrees,
        double distance_Angstrom);
} // namespace cds
#endif
