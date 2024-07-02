#ifndef INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"

#include <array>
#include <vector>

namespace cds
{
    Coordinate calculateGeometricCenter(const std::vector<Coordinate*>& coords);
    Coordinate CreateCoordinateForCenterAwayFromNeighbors(const Coordinate& centralCoord,
                                                          const std::vector<Coordinate*>& threeNeighbors,
                                                          const double distance = 1.0);
    Coordinate calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                     double angle_Degrees, double dihedral_Degrees,
                                                     double distance_Angstrom);
} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_ */
