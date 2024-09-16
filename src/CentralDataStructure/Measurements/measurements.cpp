#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CodeUtils/constants.hpp"

#include <cmath>

using cds::Coordinate;

Coordinate cds::coordinateOppositeToNeighborAverage(const Coordinate& centralCoord,
                                                    const std::vector<Coordinate>& neighbors, const double distance)
{
    Coordinate combinedVs(0.0, 0.0, 0.0);
    for (auto& neighbor : neighbors)
    {
        // normalize so that a small bond length in a H doesn't create a wonky tetrahedral
        combinedVs = combinedVs + normal(centralCoord - neighbor);
    }
    return centralCoord + scaleBy(distance, normal(combinedVs));
}

Coordinate cds::calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                      double angle_Degrees, double dihedral_Degrees,
                                                      double distanceAngstrom)
{
    double theta_Radians = constants::toRadians(angle_Degrees);
    double phi_Radians   = constants::toRadians(dihedral_Degrees);

    Coordinate lmn_y = normal(crossProduct(a - b, b - c));
    Coordinate lmn_z = normal(b - c);
    Coordinate lmn_x = crossProduct(lmn_z, lmn_y);

    double x_p = std::sin(theta_Radians) * std::cos(phi_Radians);
    double y_p = std::sin(theta_Radians) * std::sin(phi_Radians);
    double z_p = std::cos(theta_Radians);

    return c + scaleBy(distanceAngstrom, scaleBy(x_p, lmn_x) + scaleBy(y_p, lmn_y) + scaleBy(z_p, lmn_z));
}
