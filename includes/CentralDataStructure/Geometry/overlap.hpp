#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/MolecularMetadata/elements.hpp"

namespace cds
{
    int compareOverlaps(double a, double b);
    double overlapAmount(const MolecularMetadata::PotentialFactor& factor, double tolerance, const Sphere& a,
                         const Sphere& b);
    double overlapVectorSum(const std::vector<double>& vec);
    std::vector<double> overlapAboveThreshold(double threshold, const std::vector<double>& vec);
    bool containsOverlapExceedingThreshold(double threshold, const std::vector<double>& vec);
    void addOverlapsTo(std::vector<double>& vec, const std::vector<double>& added);

    inline double overlapAboveThresholdSum(double threshold, const std::vector<double>& vec)
    {
        return overlapVectorSum(overlapAboveThreshold(threshold, vec));
    }
} // namespace cds

#endif
