#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"

#include <vector>
#include <utility>

namespace cds
{
    void insertIndicesOfIntersection(std::vector<size_t>& result, double overlapTolerance, const cds::Sphere& sphere,
                                     const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    std::vector<size_t> intersectingIndices(double overlapTolerance, const cds::Sphere& sphere,
                                            const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices);
    Overlap CountOverlappingAtoms(const codeUtils::SparseVector<double>& elementRadii, double overlapTolerance,
                                  const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    void addResidueOverlaps(std::vector<Overlap>& result, const MolecularMetadata::PotentialTable& potential,
                            double overlapTolerance, const assembly::Graph& graph, const assembly::Bounds& bounds,
                            const assembly::Selection& selectionA, const assembly::Selection& selectionB,
                            const std::vector<MolecularMetadata::Element>& atomElements,
                            const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
                            size_t residueA, size_t residueB);
    std::vector<Overlap>
    overlapsBetweenSelections(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                              const assembly::Graph& graph, const assembly::Bounds& bounds,
                              const assembly::Selection& selectionA, const assembly::Selection& selectionB,
                              const std::vector<MolecularMetadata::Element>& atomElements,
                              const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge);
    std::vector<Overlap>
    overlapsWithinSelection(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                            const assembly::Graph& graph, const assembly::Bounds& bounds,
                            const assembly::Selection& selection,
                            const std::vector<MolecularMetadata::Element>& atomElements,
                            const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge);
} // namespace cds
#endif
