#ifndef INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"

#include <array>
#include <vector>

namespace cds
{
    Coordinate calculateGeometricCenter(const std::vector<Coordinate*>& coords);
    double CalculateAngle(const std::array<Coordinate*, 3>& coords, bool returnRadians = false);
    double CalculateDihedralAngle(const std::array<Coordinate*, 4>& coords, bool returnRadians = false);
    Coordinate CreateCoordinateForCenterAwayFromNeighbors(const Coordinate& centralCoord,
                                                          const std::vector<Coordinate*>& threeNeighbors,
                                                          const double distance = 1.0);
    Coordinate calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                     double angle_Degrees, double dihedral_Degrees,
                                                     double distance_Angstrom);
} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_MEASUREMENTS_MEASUREMENTS_HPP_ */
