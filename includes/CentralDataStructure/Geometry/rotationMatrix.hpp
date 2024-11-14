#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_ROTATIONMATRIX_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_ROTATIONMATRIX_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CodeUtils/references.hpp"

#include <array>
#include <vector>

namespace cds
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

        void rotateCoordinates(std::vector<CoordinateReference> coords);

      private:
        RotationMatrixInternals matrix_;
    };

    RotationMatrix rotationAroundPoint(const Coordinate& point, const Coordinate& axis, double angle);
} // namespace cds
#endif
