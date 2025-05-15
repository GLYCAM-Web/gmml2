#ifndef INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_EDITORS_SUPERIMPOSITION_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
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
    std::vector<Coordinate> matrixCoordinates(const Eigen::Matrix3Xd& matrix);
    AffineTransform affineTransform(const std::vector<Coordinate>& target, const std::vector<Coordinate>& moving);
    void Superimpose(std::vector<Coordinate>& moving, const std::vector<Coordinate>& target);
    void Superimpose(std::vector<Coordinate>& moving, const std::vector<Coordinate>& target,
                     std::vector<Coordinate>& alsoMoving);
} // namespace cds

#endif
