#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <cmath>
#include <vector>

int cds::compareOverlaps(double a, double b)
{
    return (std::fabs(a - b) <= 1e-10) ? 0 : ((a > b) ? 1 : -1);
}

double cds::overlapAmount(const MolecularMetadata::PotentialFactor& factor, double tolerance, const Sphere& a,
                          const Sphere& b)
{
    double cutoff = std::max(0.0, a.radius + b.radius - tolerance);
    double sqDist = squaredDistance(a.center, b.center);
    if (sqDist < cutoff * cutoff)
    {
        return std::max(0.0, MolecularMetadata::lennardJonesPotential(factor, sqDist));
    }
    else
    {
        return 0.0;
    }
}

double cds::overlapVectorSum(const std::vector<double>& vec)
{
    return 0.5 * codeUtils::vectorSum(0.0, vec);
}

std::vector<double> cds::overlapAboveThreshold(double threshold, const std::vector<double>& vec)
{
    std::vector<double> result;
    result.reserve(vec.size());
    for (auto& a : vec)
    {
        result.push_back(std::max(0.0, a - threshold));
    }
    return result;
}

bool cds::containsOverlapExceedingThreshold(double threshold, const std::vector<double>& vec)
{
    for (auto& a : vec)
    {
        if (a > threshold)
        {
            return true;
        }
    }
    return false;
}

void cds::addOverlapsTo(std::vector<double>& vec, const std::vector<double>& added)
{
    for (size_t n = 0; n < vec.size(); n++)
    {
        vec[n] += added[n];
    }
}
