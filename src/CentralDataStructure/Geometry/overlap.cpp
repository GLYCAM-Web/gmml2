#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"

#include <cmath>
#include <vector>

int cds::compareOverlaps(const Overlap& a, const Overlap& b)
{
    if (std::round(a.count) == std::round(b.count))
    {
        return (std::fabs(a.weight - b.weight) <= 1e-10) ? 0 : ((a.weight > b.weight) ? 1 : -1);
    }
    else
    {
        return a.count - b.count;
    }
}

cds::Overlap cds::overlapAmount(const OverlapProperties properties, const Sphere& a, const Sphere& b)
{
    double cutoff     = a.radius + b.radius - properties.tolerance;
    double sqDist     = squaredDistance(a.center, b.center);
    double weightBase = properties.weightBase;
    if (sqDist < cutoff * cutoff)
    {
        return Overlap {1.0, weightBase / (weightBase + sqDist)};
    }
    else
    {
        return Overlap {0.0, 0.0};
    }
}
