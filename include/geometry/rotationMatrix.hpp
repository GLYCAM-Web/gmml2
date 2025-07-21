#ifndef INCLUDES_GEOMETRY_ROTATIONMATRIX_HPP
#define INCLUDES_GEOMETRY_ROTATIONMATRIX_HPP

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <array>
#include <vector>

namespace gmml
{
    typedef std::array<std::array<double, 4>, 3> RotationMatrixInternals;

    class RotationMatrix
    {
      public:
        RotationMatrix(const RotationMatrixInternals& mat) : matrix_(mat) {};

        Coordinate operator*(const Coordinate& coord) const
        {
            auto col = [&](int n)
            {
                const std::array<double, 4>& m = matrix_[n];
                return dotProduct(coord, {m[0], m[1], m[2]}) + m[3];
            };
            return {col(0), col(1), col(2)};
        };

        void rotateCoordinates(std::vector<Coordinate>& coords);

      private:
        RotationMatrixInternals matrix_;
    };

    RotationMatrix rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle);
} // namespace gmml

#endif
