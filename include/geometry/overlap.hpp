#ifndef INCLUDE_GEOMETRY_OVERLAP_HPP
#define INCLUDE_GEOMETRY_OVERLAP_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"

namespace gmml
{
    int compareOverlaps(double a, double b);
    double overlapAmount(const PotentialFactor& factor, double tolerance, const Sphere& a, const Sphere& b);
    double overlapVectorSum(const std::vector<double>& vec);
    std::vector<double> overlapAboveThreshold(double threshold, const std::vector<double>& vec);
    bool containsOverlapExceedingThreshold(double threshold, const std::vector<double>& vec);
    void addOverlapsTo(std::vector<double>& vec, const std::vector<double>& added);

    inline double overlapAboveThresholdSum(double threshold, const std::vector<double>& vec)
    {
        return overlapVectorSum(overlapAboveThreshold(threshold, vec));
    }
} // namespace gmml

#endif
