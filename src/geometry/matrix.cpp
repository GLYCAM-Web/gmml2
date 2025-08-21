#include "include/geometry/matrix.hpp"

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <cmath>
#include <vector>

namespace gmml
{
    void transformCoordinates(const Matrix4x4& matrix, std::vector<Coordinate>& coords)
    {
        for (auto& coord : coords)
        {
            coord = matrix * coord;
        }
    }

    Matrix4x4 translation(const Coordinate& offset)
    {
        return Matrix4x4({
            {{1, 0, 0, offset.nth(0)}, {0, 1, 0, offset.nth(1)}, {0, 0, 1, offset.nth(2)}, {0, 0, 0, 1}}
        });
    }

    Matrix4x4 rotation(const Coordinate& axis, double angle)
    {
        double cosA = std::cos(angle);
        double sinA = std::sin(angle);

        auto diag = [&](double a) { return cosA + a * a * (1 - cosA); };
        auto cell = [&](double a, double b, double c) { return a * b * (1 - cosA) + c * sinA; };

        double x = axis.nth(0);
        double y = axis.nth(1);
        double z = axis.nth(2);

        return Matrix4x4({
            {{diag(x), cell(x, y, -z), cell(x, z, y), 0},
             {cell(x, y, z), diag(y), cell(y, z, -x), 0},
             {cell(x, z, -y), cell(y, z, x), diag(z), 0},
             {0, 0, 0, 1}}
        });
    }

    Matrix4x4 rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle)
    {
        return translation(point) * rotation(axis, angle) * translation(scaleBy(-1, point));
    }
} // namespace gmml
