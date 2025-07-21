#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containerTypes.hpp"

#include <utility>
#include <vector>

namespace gmml
{
    void insertIndicesOfIntersection(
        std::vector<size_t>& result,
        double overlapTolerance,
        const Sphere& sphere,
        const std::vector<Sphere>& coords,
        const std::vector<size_t>& indices);

    std::vector<size_t> intersectingIndices(
        double overlapTolerance,
        const Sphere& sphere,
        const std::vector<Sphere>& coords,
        const std::vector<size_t>& indices);

    double CountOverlappingAtoms(
        const util::SparseVector<double>& elementRadii,
        double overlapTolerance,
        const std::vector<Atom*>& atomsA,
        const std::vector<Atom*>& atomsB);

    void addResidueOverlaps(
        std::vector<double>& result,
        const PotentialTable& potential,
        double overlapTolerance,
        const assembly::Graph& graph,
        const assembly::Bounds& bounds,
        const assembly::Selection& selectionA,
        const assembly::Selection& selectionB,
        const std::vector<Element>& atomElements,
        const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
        size_t residueA,
        size_t residueB);

    std::vector<double> overlapsBetweenSelections(
        const PotentialTable& potential,
        double overlapTolerance,
        const assembly::Graph& graph,
        const assembly::Bounds& bounds,
        const assembly::Selection& selectionA,
        const assembly::Selection& selectionB,
        const std::vector<Element>& atomElements,
        const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge);

    std::vector<double> overlapsWithinSelection(
        const PotentialTable& potential,
        double overlapTolerance,
        const assembly::Graph& graph,
        const assembly::Bounds& bounds,
        const assembly::Selection& selection,
        const std::vector<Element>& atomElements,
        const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge);
} // namespace gmml
#endif
