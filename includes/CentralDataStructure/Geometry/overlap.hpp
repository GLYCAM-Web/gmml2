#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

namespace cds
{
    struct OverlapProperties
    {
        double weightBase;
        double tolerance;
    };

    struct Overlap
    {
        double count  = 0.0;
        double weight = 0.0;

        inline Overlap operator+(const Overlap& a) const
        {
            return {count + a.count, weight + a.weight};
        }

        inline Overlap operator*(double a) const
        {
            return {count * a, weight * a};
        }

        inline Overlap& operator+=(const Overlap& a)
        {
            *this = (*this + a);
            return *this;
        }
    };

    int compareOverlaps(const Overlap& a, const Overlap& b);
    Overlap overlapAmount(const OverlapProperties properties, const Sphere& a, const Sphere& b);
    Overlap overlapVectorSum(const std::vector<Overlap>& vec);
    std::vector<Overlap> scaledOverlaps(double scale, const std::vector<Overlap>& vec);
    void addOverlapsTo(std::vector<Overlap>& vec, const std::vector<Overlap>& added);
} // namespace cds

#endif
