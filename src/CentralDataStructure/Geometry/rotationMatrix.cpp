#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CodeUtils/references.hpp"

#include <vector>
#include <cmath>

using cds::RotationMatrix;

// Functions
void RotationMatrix::rotateCoordinates(std::vector<CoordinateReference> coords)
{
    for (auto& coord : coords)
    {
        coord.set((*this) * coord.get());
    }
}

RotationMatrix cds::rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle)
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
