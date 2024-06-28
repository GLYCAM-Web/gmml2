#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_

#include "includes/CentralDataStructure/coordinate.hpp"
#include <array>
#include <vector>

namespace cds
{
    class RotationMatrix
    {
      public:
        // Constructors.
        RotationMatrix() =
            delete; // I only want the class to be instantiated in a certain way. These deleted ctors allow this.
        RotationMatrix(const RotationMatrix&)            = delete;
        RotationMatrix& operator=(const RotationMatrix&) = delete;
        RotationMatrix(Coordinate* direction, Coordinate* parent, double angle); // the way
        // Functions
        void rotateCoordinates(std::vector<Coordinate*> coords);

      private:
        std::array<std::array<double, 4>, 3> matrix_;
    };

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ROTATIONMATRIX_HPP_ */
