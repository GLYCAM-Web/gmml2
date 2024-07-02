#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP

#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include <vector>

using cds::Coordinate;

namespace cds
{
    void Superimpose(std::vector<Coordinate*>& moving, const std::vector<Coordinate*>& target);
    void Superimpose(std::vector<Coordinate*>& moving, const std::vector<Coordinate*>& target,
                     std::vector<Coordinate*>& alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
} // namespace cds

#endif // SUPERIMPOSITION_HPP
