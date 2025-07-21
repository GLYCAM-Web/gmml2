#ifndef INCLUDES_GEOMETRY_MEASUREMENTS_HPP
#define INCLUDES_GEOMETRY_MEASUREMENTS_HPP

#include "include/geometry/geometryTypes.hpp"

#include <array>
#include <vector>

namespace gmml
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
} // namespace gmml

#endif
