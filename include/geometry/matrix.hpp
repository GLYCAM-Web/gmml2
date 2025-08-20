#ifndef INCLUDE_GEOMETRY_MATRIX_HPP
#define INCLUDE_GEOMETRY_MATRIX_HPP

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <array>
#include <vector>

namespace gmml
{
    typedef std::array<std::array<double, 4>, 4> MatrixInternals4x4;

    class Matrix4x4
    {
      public:
        Matrix4x4(const MatrixInternals4x4& mat) : values(mat) {};

        inline std::array<double, 4> operator*(const std::array<double, 4>& arr) const
        {
            auto col = [&](int n)
            {
                const std::array<double, 4>& m = values[n];
                return m[0] * arr[0] + m[1] * arr[1] + m[2] * arr[2] + m[3] * arr[3];
            };
            return {col(0), col(1), col(2), col(3)};
        };

        inline Coordinate operator*(const Coordinate& coord) const
        {
            std::array<double, 4> arr = operator*({coord.nth(0), coord.nth(1), coord.nth(2), 1});
            return Coordinate(arr[0], arr[1], arr[2]);
        };

      private:
        MatrixInternals4x4 values;
    };

    Matrix4x4 rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle);

    void transformCoordinates(const Matrix4x4& matrix, std::vector<Coordinate>& coords);
} // namespace gmml

#endif
