#include "include/geometry/matrix.hpp"

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace gmml
{
    std::vector<Coordinate> transform(const Matrix4x4& matrix, const std::vector<Coordinate>& coords)
    {
        std::function<Coordinate(const Coordinate&)> func = [&](const Coordinate& coord) { return matrix * coord; };
        return util::vectorMap(func, coords);
    }

    Matrix4x4 identityMatrix() { return translation({0, 0, 0}); }

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

    Matrix4x4 superimposition(const std::array<Coordinate, 3>& target, const std::array<Coordinate, 3>& superimposed)
    {
        auto transform = [](const Matrix4x4& m, const std::array<Coordinate, 3>& coords) {
            return std::array<Coordinate, 3> {m * coords[0], m * coords[1], m * coords[2]};
        };

        Coordinate axis(0, 0, 1);
        std::pair<Matrix4x4, Matrix4x4> identities {identityMatrix(), identityMatrix()};
        std::pair<Matrix4x4, Matrix4x4> flips {
            rotation({1, 0, 0}, constants::PI_RADIAN), rotation({1, 0, 0}, -constants::PI_RADIAN)};
        auto axialRotation = [&](const Coordinate& coord)
        {
            // we're looking for a rotation that brings coord to [0, 0, 1]
            // handle the case where their cross products are zero
            if (coord.nth(0) == 0 && coord.nth(1) == 0)
            {
                return (coord.nth(2) > 0) ? identities : flips;
            }
            else
            {
                Coordinate ax = normal(crossProduct(coord, axis));
                double a = angle(axis, coord);
                return std::pair<Matrix4x4, Matrix4x4> {rotation(ax, a), rotation(ax, -a)};
            }
        };
        Matrix4x4 centerTarget = translation(scaleBy(-1, target[0]));
        Matrix4x4 s0translation = translation(scaleBy(-1, superimposed[0]));

        auto tc = transform(centerTarget, target);

        std::pair<Matrix4x4, Matrix4x4> t1transforms = axialRotation(tc[1]);
        auto tca = transform(t1transforms.first, tc);

        auto sc = transform(s0translation, superimposed);
        std::pair<Matrix4x4, Matrix4x4> s1transforms = axialRotation(sc[1]);
        Matrix4x4 s1rotation = s1transforms.first;
        auto sca = transform(s1rotation, sc);

        auto planarAngle = [](const Coordinate& coord) { return atan2(coord.nth(1), coord.nth(0)); };
        double angle2 = planarAngle(tca[2]) - planarAngle(sca[2]);

        return translation(target[0]) * t1transforms.second * rotation(axis, angle2) * s1rotation * s0translation;
    }
} // namespace gmml
