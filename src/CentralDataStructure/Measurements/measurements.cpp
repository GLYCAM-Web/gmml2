#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CodeUtils/constants.hpp"

using cds::Coordinate;

Coordinate cds::calculateGeometricCenter(const std::vector<Coordinate*>& coords)
{
    if (coords.size() == 0)
    {
        throw std::runtime_error("Oliver what were you thinking in cds::calculateGeometricCenter?");
    }
    Coordinate center(0.0, 0.0, 0.0);
    for (auto& coord : coords)
    {
        center = center + *coord;
    }
    return scaleBy(1.0 / coords.size(), center);
}

double cds::CalculateAngle(const std::array<Coordinate*, 3>& coords, bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    Coordinate b1 = *coords[0] - *coords[1];
    Coordinate b2 = *coords[2] - *coords[1];
    double angle  = acos((dotProduct(b1, b2)) / (length(b1) * length(b2) + constants::DIST_EPSILON));
    return angle * (returnRadians ? 1.0 : (180 / constants::PI_RADIAN));
}

double cds::CalculateDihedralAngle(const std::array<Coordinate*, 4>& coords, bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    Coordinate b1 = *coords[1] - *coords[0];
    Coordinate b2 = *coords[2] - *coords[1];
    Coordinate b3 = *coords[3] - *coords[2];

    Coordinate b2xb3    = crossProduct(b2, b3);
    Coordinate b1_m_b2n = scaleBy(length(b2), b1);
    Coordinate b1xb2    = crossProduct(b1, b2);

    double angle = atan2(dotProduct(b1_m_b2n, b2xb3), dotProduct(b1xb2, b2xb3));
    return angle * (returnRadians ? 1.0 : (180 / constants::PI_RADIAN));
}

Coordinate cds::CreateCoordinateForCenterAwayFromNeighbors(const Coordinate& centralCoord,
                                                           const std::vector<Coordinate*>& threeNeighbors,
                                                           const double distance)
{
    Coordinate combinedVs(0.0, 0.0, 0.0);
    for (auto& neighbor : threeNeighbors)
    {
        // normalize so that a small bond length in a H doesn't create a wonky tetrahedral
        combinedVs = combinedVs + normal(centralCoord - *neighbor);
    }
    return centralCoord + scaleBy(distance, combinedVs);
}

Coordinate cds::calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                      double angle_Degrees, double dihedral_Degrees,
                                                      double distanceAngstrom)
{
    double theta_Radians = constants::degree2Radian(angle_Degrees);
    double phi_Radians   = constants::degree2Radian(dihedral_Degrees);

    Coordinate lmn_y = normal(crossProduct(a - b, b - c));
    Coordinate lmn_z = normal(b - c);
    Coordinate lmn_x = crossProduct(lmn_z, lmn_y);

    double x_p = sin(theta_Radians) * cos(phi_Radians);
    double y_p = sin(theta_Radians) * sin(phi_Radians);
    double z_p = cos(theta_Radians);

    return c + scaleBy(distanceAngstrom, scaleBy(x_p, lmn_x) + scaleBy(y_p, lmn_y) + scaleBy(z_p, lmn_z));
}
