#ifndef INCLUDES_CENTRALDATASTRUCTURE_ORIENTATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ORIENTATION_HPP

#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"

#include <array>

namespace cds
{
    typedef std::array<Coordinate, 4> DihedralCoordinates;

    Coordinate axis(const DihedralCoordinates& coords);
    Coordinate axis(const std::array<Coordinate, 3>& coords);
    double angle(const DihedralCoordinates& coords);
    double angle(const std::array<Coordinate, 3>& coords);
    RotationMatrix rotationTo(const DihedralCoordinates& coords, double targetAngle);
    RotationMatrix rotationTo(const std::array<Coordinate, 3>& coords, double targetAngle);
} // namespace cds

#endif
