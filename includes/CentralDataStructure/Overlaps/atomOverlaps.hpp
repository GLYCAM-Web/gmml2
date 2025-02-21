#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"

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

    void insertIndicesOfIntersection(std::vector<size_t>& result, double overlapTolerance, const cds::Sphere& sphere,
                                     const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    std::vector<size_t> intersectingIndices(double overlapTolerance, const cds::Sphere& sphere,
                                            const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    Overlap CountOverlappingAtoms(OverlapProperties properties, const std::vector<Atom*>& atomsA,
                                  const std::vector<Atom*>& atomsB);
    std::vector<Overlap>
    CountOverlappingAtoms(const MolecularMetadata::PotentialTable& potential, OverlapProperties properties,
                          const assembly::Graph& graph, const assembly::Bounds& bounds,
                          const std::vector<double>& residueWeights,
                          const std::vector<MolecularMetadata::Element>& atomElements,
                          const std::vector<bool>& includedAtoms, const std::vector<BondedResidueOverlapInput>& bonds,
                          const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB);
} // namespace cds
#endif
