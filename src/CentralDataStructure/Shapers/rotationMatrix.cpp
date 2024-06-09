#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CodeUtils/constants.hpp"

#include <cmath>

using cds::RotationMatrix;

// Functions
void RotationMatrix::rotateCoordinates(std::vector<Coordinate*> coords)
{
    for (auto& coord : coords)
    {
        *coord = (*this) * (*coord);
    }
}

RotationMatrix cds::axisAnglePositionToMatrix(Coordinate axis, double angle, const Coordinate& position)
{
    RotationMatrixInternals mat {
        {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
    };
    ;
    double u = axis.GetX();
    double v = axis.GetY();
    double w = axis.GetZ();

    double a = position.GetX();
    double b = position.GetY();
    double c = position.GetZ();

    double u2   = u * u;
    double v2   = v * v;
    double w2   = w * w;
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);

    mat[0][3] =
        a * (v2 + w2) - u * (b * v + c * w) + (u * (b * v + c * w) - a * (v2 + w2)) * cosA + (b * w - c * v) * sinA;
    mat[1][3] =
        b * (u2 + w2) - v * (a * u + c * w) + (v * (a * u + c * w) - b * (u2 + w2)) * cosA + (c * u - a * w) * sinA;
    mat[2][3] =
        c * (u2 + v2) - w * (a * u + b * v) + (w * (a * u + b * v) - c * (u2 + v2)) * cosA + (a * v - b * u) * sinA;

    mat[0][0] = u2 + (v2 + w2) * cosA;
    mat[0][1] = u * v * (1 - cosA) - w * sinA;
    mat[0][2] = u * w * (1 - cosA) + v * sinA;

    mat[1][0] = u * v * (1 - cosA) + w * sinA;
    mat[1][1] = v2 + (u2 + w2) * cosA;
    mat[1][2] = v * w * (1 - cosA) - u * sinA;

    mat[2][0] = u * w * (1 - cosA) - v * sinA;
    mat[2][1] = v * w * (1 - cosA) + u * sinA;
    mat[2][2] = w2 + (u2 + v2) * cosA;

    return RotationMatrix {mat};
}

RotationMatrix cds::dihedralToMatrix(const std::array<Coordinate*, 4> dihedral, const double dihedral_angle)
{
    double current_dihedral = CalculateDihedralAngle(dihedral, true);
    Coordinate axis         = normal(*dihedral[1] - *dihedral[2]);
    double angle            = current_dihedral - constants::degree2Radian(dihedral_angle);
    return axisAnglePositionToMatrix(axis, angle, *dihedral[1]);
}

RotationMatrix cds::angleToMatrix(const std::array<Coordinate*, 3> coords, const double angle)
{
    double current_angle  = cds::CalculateAngle(coords, true);
    double rotation_angle = constants::degree2Radian(angle) - current_angle;
    Coordinate axis       = normal(crossProduct(*coords[0] - *coords[1], *coords[2] - *coords[1]));
    return axisAnglePositionToMatrix(axis, rotation_angle, *coords[1]);
}
