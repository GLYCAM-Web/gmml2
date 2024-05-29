#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CodeUtils/constants.hpp"
#include <eigen3/Eigen/Geometry>

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

double cds::CalculateAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, const bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    double current_angle = 0.0;
    Coordinate b1        = *a1 - *a2;
    Coordinate b2        = *a3 - *a2;
    current_angle        = acos((dotProduct(b1, b2)) / (length(b1) * length(b2) + constants::DIST_EPSILON));
    if (returnRadians)
    {
        return current_angle;
    }
    return (current_angle * (180 / constants::PI_RADIAN)); // Convert to degrees
}

double cds::CalculateDihedralAngle(Coordinate* a1, Coordinate* a2, Coordinate* a3, Coordinate* a4,
                                   const bool returnRadians)
{ // returns Degrees by default, must set bool returnRadians to true for radians.
    Coordinate b1 = *a2 - *a1;
    Coordinate b2 = *a3 - *a2;
    Coordinate b3 = *a4 - *a3;
    //    Coordinate b4 = b2;
    //    b4.operator*(-1);

    Coordinate b2xb3 = crossProduct(b2, b3);

    Coordinate b1_m_b2n = scaleBy(length(b2), b1);

    Coordinate b1xb2 = crossProduct(b1, b2);

    double current_dihedral_angle = atan2(dotProduct(b1_m_b2n, b2xb3), dotProduct(b1xb2, b2xb3));
    if (returnRadians)
    {
        return current_dihedral_angle;
    }
    return (current_dihedral_angle * (180 / constants::PI_RADIAN)); // Convert to degrees
}

Coordinate cds::CreateCoordinateForCenterAwayFromNeighbors(const Coordinate* centralCoord,
                                                           std::vector<Coordinate*> threeNeighbors,
                                                           const double distance)
{
    Coordinate combinedVs(0.0, 0.0, 0.0);
    for (auto& neighbor : threeNeighbors)
    {
        combinedVs =
            combinedVs +
            normal(*centralCoord -
                   *neighbor); // normalize so that a small bond length in a H doesn't create a wonky tetrahedral
    }
    return *centralCoord + scaleBy(distance, combinedVs);
}

Coordinate cds::calculateCoordinateFromInternalCoords(const Coordinate& a, const Coordinate& b, const Coordinate& c,
                                                      double angle_Degrees, double dihedral_Degrees,
                                                      double distance_Angstrom)
{
    //  std::cout << "Distance: " << distance_Angstrom << std::endl;
    //  std::cout << "Angle: " << angle_Degrees << std::endl;
    //  std::cout << "Dihedral: " << dihedral_Degrees << std::endl;
    double theta_Radians = constants::degree2Radian(angle_Degrees);
    double phi_Radians   = constants::degree2Radian(dihedral_Degrees);

    // ToDo no. Overload the - operator properly in Coordinate.
    Coordinate cb = b - c; // original
    Coordinate ba = a - b; // original

    Coordinate lmn_y = normal(crossProduct(ba, cb));
    Coordinate lmn_z = normal(cb);
    Coordinate lmn_x = crossProduct(lmn_z, lmn_y);

    double x_p = distance_Angstrom * sin(theta_Radians) * cos(phi_Radians);
    double y_p = distance_Angstrom * sin(theta_Radians) * sin(phi_Radians);
    double z_p = distance_Angstrom * cos(theta_Radians);

    return scaleBy(x_p, lmn_x) + scaleBy(y_p, lmn_y) + scaleBy(z_p, lmn_z) + c;
}