#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CodeUtils/references.hpp"

#include <vector>

namespace cds
{
    void Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target);
    void Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target,
                     std::vector<CoordinateReference>& alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
} // namespace cds

#endif
