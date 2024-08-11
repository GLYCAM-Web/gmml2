#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>
#include <cmath>

namespace cds
{
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

    inline int compareOverlaps(const Overlap& a, const Overlap& b)
    {
        if (a.count == b.count)
        {
            return (std::fabs(a.weight - b.weight) <= 1e-10) ? 0 : ((a.weight > b.weight) ? 1 : -1);
        }
        else
        {
            return a.count - b.count;
        }
    }

    struct ResiduesWithOverlapWeight
    {
        std::vector<Residue*> residues;
        std::vector<double> weights;
    };

    struct ResidueAtomOverlapInput
    {
        std::vector<Sphere> atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<std::pair<size_t, size_t>> residueAtoms;
        const std::vector<double> residueWeights;
    };

    struct ResidueAtomOverlapInputReference
    {
        std::vector<Sphere>& atomCoordinates;
        std::vector<Sphere>& boundingSpheres;
        const std::vector<std::pair<size_t, size_t>>& residueAtoms;
        const std::vector<double>& residueWeights;
    };

    ResidueAtomOverlapInput toOverlapInput(const ResiduesWithOverlapWeight& input);
    Overlap CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    Overlap CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                  const ResidueAtomOverlapInputReference& moving);
    Overlap CountOverlappingAtoms(const ResiduesWithOverlapWeight& residuesA,
                                  const ResiduesWithOverlapWeight& residuesB);
} // namespace cds
#endif
