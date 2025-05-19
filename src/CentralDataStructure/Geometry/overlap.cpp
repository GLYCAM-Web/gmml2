#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/MolecularMetadata/elements.hpp"

#include <cmath>
#include <vector>

int cds::compareOverlaps(const Overlap& a, const Overlap& b)
{
    return (std::fabs(a.weight - b.weight) <= 1e-10) ? 0 : ((a.weight > b.weight) ? 1 : -1);
}

cds::Overlap cds::overlapAmount(const MolecularMetadata::PotentialFactor& factor, double tolerance, const Sphere& a,
                                const Sphere& b)
{
    double cutoff = std::max(0.0, a.radius + b.radius - tolerance);
    double sqDist = squaredDistance(a.center, b.center);
    if (sqDist < cutoff * cutoff)
    {
        return Overlap {1.0, MolecularMetadata::lennardJonesPotential(factor, sqDist)};
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

bool cds::containsOverlapExceedingThreshold(double threshold, const std::vector<cds::Overlap>& vec)
{
    for (auto& a : vec)
    {
        if (a.weight > threshold)
        {
            return true;
        }
    }
    return false;
}

void cds::scaleOverlaps(double scale, std::vector<Overlap>& vec)
{
    for (size_t n = 0; n < vec.size(); n++)
    {
        vec[n] = vec[n] * scale;
    }
}

std::vector<cds::Overlap> cds::scaledOverlaps(double scale, const std::vector<Overlap>& vec)
{
    std::vector<cds::Overlap> result = vec;
    scaleOverlaps(scale, result);
    return result;
}

void cds::addOverlapsTo(std::vector<Overlap>& vec, const std::vector<Overlap>& added)
{
    for (size_t n = 0; n < vec.size(); n++)
    {
        vec[n] += added[n];
    }
}
