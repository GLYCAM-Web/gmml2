#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
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

    struct BondedResidueOverlapInput
    {
        std::array<size_t, 2> residueIndices;
        std::array<std::vector<bool>, 2> ignoredAtoms;
    };

    void insertIndicesOfIntersection(std::vector<size_t>& result, const cds::Sphere& sphere,
                                     const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    std::vector<size_t> intersectingIndices(const cds::Sphere& sphere, const std::vector<cds::Sphere>& coords,
                                            const std::vector<size_t>& indices);
    Overlap CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    Overlap CountOverlappingAtoms(const std::vector<Sphere>& atomBounds, const std::vector<Sphere>& residueBounds,
                                  const std::vector<std::vector<size_t>>& residueAtoms,
                                  const std::vector<double>& residueWeights, const std::vector<bool>& ignoredAtoms,
                                  const std::vector<BondedResidueOverlapInput>& bonds,
                                  const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB);
} // namespace cds
#endif
