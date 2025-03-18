#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/references.hpp"
#include <eigen3/Eigen/Geometry>

#include <vector>

namespace cds
{
    struct AffineTransform
    {
        Eigen::Matrix3Xd moving;
        Eigen::Matrix3Xd target;
        Eigen::Affine3d affine;
    };

    Eigen::Matrix3Xd generateMatrix(const std::vector<Coordinate>& coordinates);
    Eigen::Matrix3Xd generateMatrix(const std::vector<CoordinateReference>& coordinates);
    std::vector<Coordinate> matrixCoordinates(const Eigen::Matrix3Xd& matrix);
    AffineTransform affineTransform(const std::vector<Coordinate>& target, const std::vector<Coordinate>& moving);
    void Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target);
    void Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target,
                     std::vector<CoordinateReference>& alsoMoving);
    // A function to test Find3DAffineTransform()
    // void TestFind3DAffineTransform();
} // namespace cds

#endif
