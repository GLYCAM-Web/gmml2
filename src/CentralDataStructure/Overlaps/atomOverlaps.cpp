#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"

#include <numeric>
#include <stdexcept>

namespace
{
    void insertNonIgnored(std::vector<size_t>& result, const std::vector<size_t>& indices,
                          const std::vector<bool>& includeAtom, const std::vector<bool>& localIgnore)
    {
        if (localIgnore.size() != indices.size())
        {
            throw std::runtime_error("panic");
        }
        result.reserve(indices.size());
        for (size_t n = 0; n < indices.size(); n++)
        {
            if (includeAtom[indices[n]] && !localIgnore[n])
            {
                result.push_back(indices[n]);
            }
        }
    }

    void insertIntersection(std::vector<size_t>& result, double overlapTolerance, const cds::Sphere& sphere,
                            const std::vector<cds::Sphere>& bounds, const std::vector<bool>& includeAtom,
                            const std::vector<size_t>& indices)
    {
        result.reserve(indices.size());
        for (size_t index : indices)
        {
            if (includeAtom[index] && cds::spheresOverlap(overlapTolerance, sphere, bounds[index]))
            {
                result.push_back(index);
            }
        }
    }

    std::array<std::vector<size_t>, 2>
    atomsToCheck(const assembly::Graph& graph, const assembly::Bounds& bounds, const assembly::Selection& selectionA,
                 const assembly::Selection& selectionB,
                 const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge, double overlapTolerance,
                 size_t residueA, size_t residueB)
    {
        const cds::Sphere& residueBoundsA = bounds.residues[residueA];
        const cds::Sphere& residueBoundsB = bounds.residues[residueB];
        const std::vector<size_t>& atomsA = residueAtoms(graph, residueA);
        const std::vector<size_t>& atomsB = residueAtoms(graph, residueB);
        std::vector<size_t> indicesA;
        std::vector<size_t> indicesB;
        size_t adjacency                 = codeUtils::indexOf(graph.residues.nodes.nodeAdjacencies[residueA], residueB);
        const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residueA];
        if (adjacency < edges.size())
        {
            size_t edgeIndex                                     = edges[adjacency];
            const std::array<std::vector<bool>, 2>& ignoredAtoms = residueAtomsCloseToEdge[edgeIndex];
            bool order = !(graph.residues.edges.nodeAdjacencies[edgeIndex][0] == residueA);
            insertNonIgnored(indicesA, atomsA, selectionA.atoms, ignoredAtoms[order]);
            insertNonIgnored(indicesB, atomsB, selectionB.atoms, ignoredAtoms[!order]);
        }
        else if (cds::spheresOverlap(overlapTolerance, residueBoundsA, residueBoundsB))
        {
            insertIntersection(indicesA, overlapTolerance, residueBoundsB, bounds.atoms, selectionA.atoms, atomsA);
            insertIntersection(indicesB, overlapTolerance, residueBoundsA, bounds.atoms, selectionB.atoms, atomsB);
        }
        return {indicesA, indicesB};
    }
} // namespace

void cds::insertIndicesOfIntersection(std::vector<size_t>& result, double overlapTolerance, const Sphere& sphere,
                                      const std::vector<Sphere>& coords, const std::vector<size_t>& indices)
{
    result.reserve(indices.size());
    for (size_t index : indices)
    {
        auto& a = coords[index];
        if (cds::spheresOverlap(overlapTolerance, sphere, a))
        {
            result.push_back(index);
        }
    }
}

std::vector<size_t> cds::intersectingIndices(double overlapTolerance, const cds::Sphere& sphere,
                                             const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices)
{
    std::vector<size_t> result;
    insertIndicesOfIntersection(result, overlapTolerance, sphere, coords, indices);
    return result;
}

cds::Overlap cds::CountOverlappingAtoms(const codeUtils::SparseVector<double>& elementRadii, double overlapTolerance,
                                        const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    std::vector<Sphere> coordsA                       = atomCoordinatesWithRadii(elementRadii, atomsA);
    std::vector<Sphere> coordsB                       = atomCoordinatesWithRadii(elementRadii, atomsB);
    std::vector<MolecularMetadata::Element> elementsA = cds::atomElements(atomsA);
    std::vector<MolecularMetadata::Element> elementsB = cds::atomElements(atomsB);

    Overlap overlap {0, 0.0};
    for (size_t n = 0; n < atomsA.size(); n++)
    {
        for (size_t k = 0; k < atomsB.size(); k++)
        {
            double scale =
                MolecularMetadata::potentialWeight(MolecularMetadata::potentialTable(), elementsA[n], elementsB[k]);
            overlap += overlapAmount(overlapTolerance, scale, coordsA[n], coordsB[k]);
        }
    }
    return overlap;
}

void cds::addResidueOverlaps(std::vector<Overlap>& result, const MolecularMetadata::PotentialTable& potential,
                             double overlapTolerance, const assembly::Graph& graph, const assembly::Bounds& bounds,
                             const assembly::Selection& selectionA, const assembly::Selection& selectionB,
                             const std::vector<MolecularMetadata::Element>& atomElements,
                             const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
                             size_t residueA, size_t residueB)
{
    std::array<std::vector<size_t>, 2> toCheck = atomsToCheck(
        graph, bounds, selectionA, selectionB, residueAtomsCloseToEdge, overlapTolerance, residueA, residueB);
    for (size_t atomA : toCheck[0])
    {
        MolecularMetadata::Element elementA = atomElements[atomA];
        for (size_t atomB : toCheck[1])
        {
            MolecularMetadata::Element elementB = atomElements[atomB];
            double scale                        = MolecularMetadata::potentialWeight(potential, elementA, elementB);
            Overlap overlap = overlapAmount(overlapTolerance, scale, bounds.atoms[atomA], bounds.atoms[atomB]);
            result[atomA]   += overlap;
            result[atomB]   += overlap;
        }
    }
}

std::vector<cds::Overlap>
cds::overlapsBetweenSelections(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                               const assembly::Graph& graph, const assembly::Bounds& bounds,
                               const assembly::Selection& selectionA, const assembly::Selection& selectionB,
                               const std::vector<MolecularMetadata::Element>& atomElements,
                               const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge)
{
    std::vector<Overlap> result(graph.indices.atomCount, {0.0, 0.0});
    std::vector<size_t> moleculesA                     = selectedMolecules(selectionA);
    std::vector<size_t> moleculesB                     = selectedMolecules(selectionB);
    std::vector<std::vector<size_t>> moleculeResiduesA = moleculeSelectedResidues(graph, selectionA, moleculesA);
    std::vector<std::vector<size_t>> moleculeResiduesB = moleculeSelectedResidues(graph, selectionB, moleculesB);
    for (size_t maIndex = 0; maIndex < moleculesA.size(); maIndex++)
    {
        size_t moleculeA    = moleculesA[maIndex];
        cds::Sphere boundsA = bounds.molecules[moleculeA];
        for (size_t mbIndex = 0; mbIndex < moleculesB.size(); mbIndex++)
        {
            size_t moleculeB    = moleculesB[mbIndex];
            cds::Sphere boundsB = bounds.molecules[moleculeB];
            if (cds::spheresOverlap(overlapTolerance, bounds.molecules[moleculeA], bounds.molecules[moleculeB]))
            {
                std::vector<size_t> residuesA =
                    intersectingIndices(overlapTolerance, boundsB, bounds.residues, moleculeResiduesA[maIndex]);
                std::vector<size_t> residuesB =
                    intersectingIndices(overlapTolerance, boundsA, bounds.residues, moleculeResiduesB[mbIndex]);
                for (size_t residueA : residuesA)
                {
                    for (size_t residueB : residuesB)
                    {
                        if (residueA != residueB &&
                            cds::spheresOverlap(overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                        {
                            addResidueOverlaps(result, potential, overlapTolerance, graph, bounds, selectionA,
                                               selectionB, atomElements, residueAtomsCloseToEdge, residueA, residueB);
                        }
                    }
                }
            }
        }
    }
    return result;
}

std::vector<cds::Overlap>
cds::overlapsWithinSelection(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                             const assembly::Graph& graph, const assembly::Bounds& bounds,
                             const assembly::Selection& selection,
                             const std::vector<MolecularMetadata::Element>& atomElements,
                             const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge)
{
    std::vector<Overlap> result(graph.indices.atomCount, {0.0, 0.0});
    std::vector<size_t> molecules                     = selectedMolecules(selection);
    std::vector<std::vector<size_t>> moleculeResidues = moleculeSelectedResidues(graph, selection, molecules);
    for (size_t maIndex = 0; maIndex < molecules.size(); maIndex++)
    {
        size_t moleculeA    = molecules[maIndex];
        cds::Sphere boundsA = bounds.molecules[moleculeA];
        for (size_t mbIndex = maIndex; mbIndex < molecules.size(); mbIndex++)
        {
            size_t moleculeB = molecules[mbIndex];
            if (moleculeA == moleculeB)
            {
                const std::vector<size_t>& residues = moleculeResidues[maIndex];
                for (size_t ra = 0; ra < residues.size(); ra++)
                {
                    size_t residueA = residues[ra];
                    for (size_t rb = ra + 1; rb < residues.size(); rb++)
                    {
                        size_t residueB = residues[rb];
                        if (cds::spheresOverlap(overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                        {
                            addResidueOverlaps(result, potential, overlapTolerance, graph, bounds, selection, selection,
                                               atomElements, residueAtomsCloseToEdge, residueA, residueB);
                        }
                    }
                }
            }
            else
            {
                cds::Sphere boundsB = bounds.molecules[moleculeB];
                std::vector<size_t> residuesA =
                    intersectingIndices(overlapTolerance, boundsB, bounds.residues, moleculeResidues[maIndex]);
                std::vector<size_t> residuesB =
                    intersectingIndices(overlapTolerance, boundsA, bounds.residues, moleculeResidues[mbIndex]);
                for (size_t residueA : residuesA)
                {
                    for (size_t residueB : residuesB)
                    {
                        if (cds::spheresOverlap(overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                        {
                            addResidueOverlaps(result, potential, overlapTolerance, graph, bounds, selection, selection,
                                               atomElements, residueAtomsCloseToEdge, residueA, residueB);
                        }
                    }
                }
            }
        }
    }
    return result;
}
