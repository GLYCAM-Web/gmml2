#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"
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
                auto& m = matrix_[n];
                return dotProduct(coord, {m[0], m[1], m[2]}) + m[3];
            };
            return {col(0), col(1), col(2)};
        };

        void rotateCoordinates(std::vector<Coordinate*> coords);

      private:
        RotationMatrixInternals matrix_;
    };

    RotationMatrix rotationAroundPoint(const Coordinate point, const Coordinate axis, double angle);
    RotationMatrix dihedralToMatrix(const std::array<Coordinate*, 4> dihedral, const double dihedral_angle);
    RotationMatrix angleToMatrix(const std::array<Coordinate*, 3> coords, const double angle);

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_ */
