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
    struct AtomOverlapData
    {
        const std::vector<Sphere>& bounds;
        const std::vector<MolecularMetadata::Element>& elements;
        const std::vector<bool>& included;
    };

    struct ResidueOverlapData
    {
        const std::vector<Sphere>& bounds;
        const std::vector<double>& weights;
    };

    struct ResiduesWithOverlapWeight
    {
        std::vector<Residue*> residues;
        std::vector<double> weights;
    };

    void insertIndicesOfIntersection(std::vector<size_t>& result, double overlapTolerance, const cds::Sphere& sphere,
                                     const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    std::vector<size_t> intersectingIndices(double overlapTolerance, const cds::Sphere& sphere,
                                            const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    Overlap CountOverlappingAtoms(double overlapTolerance, const std::vector<Atom*>& atomsA,
                                  const std::vector<Atom*>& atomsB);
    std::vector<Overlap>
    CountOverlappingAtoms(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                          const assembly::Graph& graph, const AtomOverlapData& atomData,
                          const ResidueOverlapData& residueData,
                          const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
                          const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB);
} // namespace cds
#endif
