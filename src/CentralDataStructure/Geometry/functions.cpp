#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"

#include <cmath>
#include <vector>

using cds::Coordinate;

double cds::length(const Coordinate& a)
{
    return std::sqrt(squaredLength(a));
}

double cds::distance(const Coordinate& a, const Coordinate& b)
{
    return std::sqrt(squaredDistance(a, b));
}

Coordinate cds::scaleBy(double factor, const Coordinate& a)
{
    auto scaled = [&](int n)
    {
        return factor * a.nth(n);
    };
    return {scaled(0), scaled(1), scaled(2)};
}

Coordinate cds::normal(const Coordinate& a)
{
    return scaleBy(1.0 / length(a), a);
}

Coordinate cds::crossProduct(const Coordinate& a, const Coordinate& b)
{
    auto cross = [&](int n, int k)
    {
        return a.nth(n) * b.nth(k) - a.nth(k) * b.nth(n);
    };
    return {cross(1, 2), cross(2, 0), cross(0, 1)};
}

Coordinate cds::projection(const Coordinate& a, const Coordinate& b)
{
    return scaleBy(dotProduct(a, b) / dotProduct(b, b), b);
}

Coordinate cds::coordinateMean(const std::vector<Coordinate>& coords)
{
    Coordinate center(0.0, 0.0, 0.0);
    for (auto& coord : coords)
    {
        center = center + coord;
    }
    return scaleBy(1.0 / coords.size(), center);
}
