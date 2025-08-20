#include "include/geometry/rotationMatrix.hpp"

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <cmath>
#include <vector>

namespace gmml
{
    void RotationMatrix::rotateCoordinates(std::vector<Coordinate>& coords)
    {
        for (auto& coord : coords)
        {
            coord = (*this) * coord;
        }
    }

    RotationMatrix rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle)
    {
        double cosA = std::cos(angle);
        double sinA = std::sin(angle);

        auto diag = [&](double a) { return cosA + a * a * (1 - cosA); };
        auto cell = [&](double a, double b, double c) { return a * b * (1 - cosA) + c * sinA; };

        double x = axis.nth(0);
        double y = axis.nth(1);
        double z = axis.nth(2);

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
} // namespace gmml
