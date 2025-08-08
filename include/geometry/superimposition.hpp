#ifndef INCLUDE_GEOMETRY_SUPERIMPOSITION_HPP
#define INCLUDE_GEOMETRY_SUPERIMPOSITION_HPP

#include "include/geometry/geometryTypes.hpp"

#include <array>
#include <eigen3/Eigen/Geometry>
#include <vector>

namespace gmml
{
    Eigen::Matrix3Xd generateMatrix(const std::vector<Coordinate>& coordinates);
    std::vector<Coordinate> matrixCoordinates(const Eigen::Matrix3Xd& matrix);
    Eigen::Affine3d affineTransform(const Eigen::Matrix3Xd& target, const Eigen::Matrix3Xd& initial);
} // namespace gmml

#endif
