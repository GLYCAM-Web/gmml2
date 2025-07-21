#include "include/geometry/orientation.hpp"

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/rotationMatrix.hpp"

#include <cmath>
#include <limits>

namespace gmml
{
    Coordinate axis(const DihedralCoordinates& coords) { return normal(coords[1] - coords[2]); }

    Coordinate axis(const std::array<Coordinate, 3>& coords)
    {
        return normal(crossProduct(coords[0] - coords[1], coords[1] - coords[2]));
    }

    double angle(const DihedralCoordinates& coords)
    {
        Coordinate b1 = coords[1] - coords[0];
        Coordinate b2 = coords[2] - coords[1];
        Coordinate b3 = coords[3] - coords[2];

        Coordinate b2xb3 = crossProduct(b2, b3);
        Coordinate b1xb2 = crossProduct(b1, b2);

        return std::atan2(dotProduct(scaleBy(length(b2), b1), b2xb3), dotProduct(b1xb2, b2xb3));
    }

    double angle(const std::array<Coordinate, 3>& coords)
    {
        Coordinate b1 = coords[0] - coords[1];
        Coordinate b2 = coords[2] - coords[1];
        return std::acos((dotProduct(b1, b2)) / (length(b1) * length(b2) + std::numeric_limits<double>::epsilon()));
    }

    RotationMatrix rotationTo(const DihedralCoordinates& coords, double targetAngle)
    {
        return rotationAroundPoint(coords[1], axis(coords), angle(coords) - targetAngle);
    }

    RotationMatrix rotationTo(const std::array<Coordinate, 3>& coords, double targetAngle)
    {
        return rotationAroundPoint(coords[1], axis(coords), angle(coords) - targetAngle);
    }
} // namespace gmml
