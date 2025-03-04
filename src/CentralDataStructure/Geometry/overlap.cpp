#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"

#include <cmath>
#include <vector>

int cds::compareOverlaps(const Overlap& a, const Overlap& b)
{
    return (std::fabs(a.weight - b.weight) <= 1e-10) ? 0 : ((a.weight > b.weight) ? 1 : -1);
}

cds::Overlap cds::overlapAmount(double tolerance, double scale, const Sphere& a, const Sphere& b)
{
    double cutoff = std::max(0.0, a.radius + b.radius - tolerance);
    double sqDist = squaredDistance(a.center, b.center);
    if (sqDist < cutoff * cutoff)
    {
        double pow4  = sqDist * sqDist;
        double pow12 = pow4 * pow4 * pow4;
        return Overlap {1.0, scale / (std::numeric_limits<double>::epsilon() + pow12)};
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

bool cds::containsOverlap(const std::vector<cds::Overlap>& vec)
{
    for (auto& a : vec)
    {
        if (a.count > 0)
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
