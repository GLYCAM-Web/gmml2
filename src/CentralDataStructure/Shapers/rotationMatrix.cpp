#include "includes/CentralDataStructure/Shapers/rotationMatrix.hpp"

#include <cmath>

using cds::RotationMatrix;

// Constructors:
RotationMatrix::RotationMatrix(Coordinate direction, const Coordinate& parent, double angle) // the way
{
    matrix_ = std::array<std::array<double, 4>, 3> {
        {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
    };
    direction = normal(direction);
    double u  = direction.GetX();
    double v  = direction.GetY();
    double w  = direction.GetZ();

    double a = parent.GetX();
    double b = parent.GetY();
    double c = parent.GetZ();

    double u2                 = u * u;
    double v2                 = v * v;
    double w2                 = w * w;
    double cos_rotation_angle = std::cos(angle);
    double sin_rotation_angle = std::sin(angle);

    matrix_[0][3] = a * (v2 + w2) - u * (b * v + c * w) + (u * (b * v + c * w) - a * (v2 + w2)) * cos_rotation_angle +
                    (b * w - c * v) * sin_rotation_angle;
    matrix_[1][3] = b * (u2 + w2) - v * (a * u + c * w) + (v * (a * u + c * w) - b * (u2 + w2)) * cos_rotation_angle +
                    (c * u - a * w) * sin_rotation_angle;
    matrix_[2][3] = c * (u2 + v2) - w * (a * u + b * v) + (w * (a * u + b * v) - c * (u2 + v2)) * cos_rotation_angle +
                    (a * v - b * u) * sin_rotation_angle;

    matrix_[0][0] = u2 + (v2 + w2) * cos_rotation_angle;
    matrix_[0][1] = u * v * (1 - cos_rotation_angle) - w * sin_rotation_angle;
    matrix_[0][2] = u * w * (1 - cos_rotation_angle) + v * sin_rotation_angle;

    matrix_[1][0] = u * v * (1 - cos_rotation_angle) + w * sin_rotation_angle;
    matrix_[1][1] = v2 + (u2 + w2) * cos_rotation_angle;
    matrix_[1][2] = v * w * (1 - cos_rotation_angle) - u * sin_rotation_angle;

    matrix_[2][0] = u * w * (1 - cos_rotation_angle) - v * sin_rotation_angle;
    matrix_[2][1] = v * w * (1 - cos_rotation_angle) + u * sin_rotation_angle;
    matrix_[2][2] = w2 + (u2 + v2) * cos_rotation_angle;
}

// Functions
void RotationMatrix::rotateCoordinates(std::vector<Coordinate*> coords)
{
    auto col = [&](const Coordinate& coord, int n)
    {
        auto& m = matrix_[n];
        return dotProduct(coord, {m[0], m[1], m[2]}) + m[3];
    };
    for (auto& coord : coords)
    {
        coord->set({col(*coord, 0), col(*coord, 1), col(*coord, 2)});
    }
}
