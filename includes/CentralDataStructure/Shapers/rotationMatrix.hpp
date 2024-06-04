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
        // Functions
        void rotateCoordinates(std::vector<Coordinate*> coords);

      private:
        RotationMatrixInternals matrix_;
    };

    RotationMatrix axisAnglePositionToMatrix(Coordinate axis, double angle, const Coordinate& position);
    RotationMatrix dihedralToMatrix(const std::array<Coordinate*, 4> dihedral, const double dihedral_angle);
    RotationMatrix angleToMatrix(const std::array<Coordinate*, 3> coords, const double angle);

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_ */
