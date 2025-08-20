#ifndef INCLUDE_GEOMETRY_ORIENTATION_HPP
#define INCLUDE_GEOMETRY_ORIENTATION_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/matrix.hpp"

#include <array>

namespace gmml
{
    typedef std::array<Coordinate, 4> DihedralCoordinates;

    Coordinate axis(const DihedralCoordinates& coords);
    Coordinate axis(const std::array<Coordinate, 3>& coords);
    double angle(const DihedralCoordinates& coords);
    double angle(const std::array<Coordinate, 3>& coords);
    Matrix4x4 rotationTo(const DihedralCoordinates& coords, double targetAngle);
    Matrix4x4 rotationTo(const std::array<Coordinate, 3>& coords, double targetAngle);
} // namespace gmml

#endif
