#ifndef INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_GEOMETRY_OVERLAP_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"

namespace cds
{
    struct OverlapProperties
    {
        double weightBase;
        double tolerance;
    };

    struct Overlap
    {
        double count;
        double weight;

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
    Overlap overlapAmount(const OverlapProperties properties, const std::vector<Sphere>& avec,
                          const std::vector<Sphere>& bvec);
} // namespace cds

#endif