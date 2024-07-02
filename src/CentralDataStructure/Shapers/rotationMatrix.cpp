#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

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

RotationMatrix cds::rotationAroundPoint(const Coordinate point, const Coordinate axis, double angle)
{
    double cosA = std::cos(angle);
    double sinA = std::sin(angle);

    auto diag = [&](double a)
    {
        return cosA + a * a * (1 - cosA);
    };
    auto cell = [&](double a, double b, double c)
    {
        return a * b * (1 - cosA) + c * sinA;
    };

    double x = axis.GetX();
    double y = axis.GetY();
    double z = axis.GetZ();

    RotationMatrixInternals mat {
        {{diag(x), cell(x, y, -z), cell(x, z, y), 0.0},
         {cell(x, y, z), diag(y), cell(y, z, -x), 0.0},
         {cell(x, z, -y), cell(y, z, x), diag(z), 0.0}}
    };

    for (size_t n = 0; n < 3; n++)
    {
        mat[n][3] = point.nth(n) - dotProduct(point, {mat[n][0], mat[n][1], mat[n][2]});
    }

    return RotationMatrix {mat};
}

RotationMatrix cds::dihedralToMatrix(const std::array<Coordinate*, 4> dihedral, const double dihedral_angle)
{
    double current_dihedral = CalculateDihedralAngle(dihedral);
    Coordinate axis         = normal(*dihedral[1] - *dihedral[2]);
    double angle            = current_dihedral - dihedral_angle;
    return rotationAroundPoint(*dihedral[1], axis, angle);
}

RotationMatrix cds::angleToMatrix(const std::array<Coordinate*, 3> coords, const double angle)
{
    double current_angle  = cds::CalculateAngle(coords);
    double rotation_angle = angle - current_angle;
    Coordinate axis       = normal(crossProduct(*coords[0] - *coords[1], *coords[2] - *coords[1]));
    return rotationAroundPoint(*coords[1], axis, rotation_angle);
}
