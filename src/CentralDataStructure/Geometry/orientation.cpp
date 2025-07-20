#include "includes/CentralDataStructure/Geometry/orientation.hpp"

#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"

#include <cmath>
#include <limits>

cds::Coordinate cds::axis(const DihedralCoordinates& coords) { return normal(coords[1] - coords[2]); }

cds::Coordinate cds::axis(const std::array<Coordinate, 3>& coords)
{
    return normal(crossProduct(coords[0] - coords[1], coords[1] - coords[2]));
}

double cds::angle(const DihedralCoordinates& coords)
{
    Coordinate b1 = coords[1] - coords[0];
    Coordinate b2 = coords[2] - coords[1];
    Coordinate b3 = coords[3] - coords[2];

    Coordinate b2xb3 = crossProduct(b2, b3);
    Coordinate b1xb2 = crossProduct(b1, b2);

    return std::atan2(dotProduct(scaleBy(length(b2), b1), b2xb3), dotProduct(b1xb2, b2xb3));
}

double cds::angle(const std::array<Coordinate, 3>& coords)
{
    Coordinate b1 = coords[0] - coords[1];
    Coordinate b2 = coords[2] - coords[1];
    return std::acos((dotProduct(b1, b2)) / (length(b1) * length(b2) + std::numeric_limits<double>::epsilon()));
}

cds::RotationMatrix cds::rotationTo(const DihedralCoordinates& coords, double targetAngle)
{
    return rotationAroundPoint(coords[1], axis(coords), angle(coords) - targetAngle);
}

cds::RotationMatrix cds::rotationTo(const std::array<Coordinate, 3>& coords, double targetAngle)
{
    return rotationAroundPoint(coords[1], axis(coords), angle(coords) - targetAngle);
}
