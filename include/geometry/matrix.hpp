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

        inline std::array<double, 4> row(int n) const { return values[n]; }

        inline std::array<double, 4> column(int n) const
        {
            return {values[0][n], values[1][n], values[2][n], values[3][n]};
        }

        inline Matrix4x4 operator*(const Matrix4x4& mat) const
        {
            MatrixInternals4x4 internals;

            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    std::array<double, 4> r = row(j);
                    std::array<double, 4> c = mat.column(k);
                    auto m = [&](int n) { return r[n] * c[n]; };
                    internals[j][k] = m(0) + m(1) + m(2) + m(3);
                }
            }

            return Matrix4x4(internals);
        }

        inline std::array<double, 4> operator*(const std::array<double, 4>& arr) const
        {
            auto col = [&](int k)
            {
                const std::array<double, 4>& row = values[k];
                auto m = [&](int n) { return row[n] * arr[n]; };
                return m(0) + m(1) + m(2) + m(3);
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

    Matrix4x4 translation(const Coordinate& offset);
    Matrix4x4 rotation(const Coordinate& axis, double angle);
    Matrix4x4 rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle);

    void transformCoordinates(const Matrix4x4& matrix, std::vector<Coordinate>& coords);
} // namespace gmml

#endif
