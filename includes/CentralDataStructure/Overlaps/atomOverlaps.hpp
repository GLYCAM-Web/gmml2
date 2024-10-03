#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    struct ResiduesWithOverlapWeight
    {
        std::vector<Residue*> residues;
        std::vector<double> weights;
    };

    struct ResidueAtomOverlapInput
    {
        std::vector<Sphere> atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<size_t> residueIndices;
        const std::vector<std::vector<size_t>> residueAtoms;
        const std::vector<double> residueWeights;
        const std::vector<bool> firstResidueBondedAtoms;
    };

    struct BondedResidueOverlapInput
    {
        std::array<size_t, 2> residueIndices;
        std::array<std::vector<bool>, 2> ignoredAtoms;
    };

    struct ResidueAtomOverlapInputReference
    {
        const std::vector<Sphere>& atomCoordinates;
        const std::vector<Sphere>& boundingSpheres;
        const std::vector<std::vector<size_t>>& residueAtoms;
        const std::vector<double>& residueWeights;
    };

    void insertIndicesOfIntersection(std::vector<size_t>& result, cds::Sphere sphere,
                                     const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    std::vector<size_t> intersectingIndices(cds::Sphere sphere, const std::vector<cds::Sphere>& coords,
                                            const std::vector<size_t>& indices);
    Overlap CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    Overlap CountOverlappingAtoms(const ResidueAtomOverlapInputReference& input,
                                  const std::vector<BondedResidueOverlapInput>& bonds,
                                  const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB);
    Overlap CountOverlappingAtoms(const ResiduesWithOverlapWeight& residuesA,
                                  const ResiduesWithOverlapWeight& residuesB);
} // namespace cds
#endif
