#include "include/geometry/overlap.hpp"

#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containers.hpp"

#include <cmath>
#include <vector>

namespace gmml
{
    int compareOverlaps(double a, double b) { return (std::fabs(a - b) <= 1e-10) ? 0 : ((a > b) ? 1 : -1); }

    double overlapAmount(const PotentialFactor& factor, double tolerance, const Sphere& a, const Sphere& b)
    {
        double cutoff = std::max(0.0, a.radius + b.radius - tolerance);
        double sqDist = squaredDistance(a.center, b.center);
        if (sqDist < cutoff * cutoff)
        {
            return std::max(0.0, lennardJonesPotential(factor, sqDist));
        }
        else
        {
            return 0.0;
        }
    }

    double overlapVectorSum(const std::vector<double>& vec) { return 0.5 * util::vectorSum(0.0, vec); }

    std::vector<double> overlapAboveThreshold(double threshold, const std::vector<double>& vec)
    {
        std::vector<double> result;
        result.reserve(vec.size());
        for (auto& a : vec)
        {
            result.push_back(std::max(0.0, a - threshold));
        }
        return result;
    }

    bool containsOverlapExceedingThreshold(double threshold, const std::vector<double>& vec)
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

    void addOverlapsTo(std::vector<double>& vec, const std::vector<double>& added)
    {
        for (size_t n = 0; n < vec.size(); n++)
        {
            vec[n] += added[n];
        }
    }
} // namespace gmml
