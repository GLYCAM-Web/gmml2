#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"

#include <cmath>
#include <vector>

int cds::compareOverlaps(const Overlap& a, const Overlap& b)
{
    return (std::fabs(a.weight - b.weight) <= 1e-10) ? 0 : ((a.weight > b.weight) ? 1 : -1);
}

cds::Overlap cds::overlapAmount(const OverlapProperties properties, const Sphere& a, const Sphere& b)
{
    double cutoff     = std::max(0.0, a.radius + b.radius - properties.tolerance);
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

cds::Overlap cds::overlapVectorSum(const std::vector<Overlap>& vec)
{
    Overlap result {0.0, 0.0};
    for (auto& a : vec)
    {
        result += a;
    }
    return result * 0.5;
}

std::vector<cds::Overlap> cds::scaledOverlaps(double scale, const std::vector<Overlap>& vec)
{
    std::vector<cds::Overlap> result;
    result.reserve(vec.size());
    for (size_t n = 0; n < vec.size(); n++)
    {
        result[n] = vec[n] * scale;
    }
    return result;
}

void cds::addOverlapsTo(std::vector<Overlap>& vec, const std::vector<Overlap>& added)
{
    for (size_t n = 0; n < vec.size(); n++)
    {
        vec[n] += added[n];
    }
}
