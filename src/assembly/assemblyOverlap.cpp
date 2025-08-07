#include "include/assembly/assemblyOverlap.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containers.hpp"

#include <numeric>
#include <stdexcept>
#include <vector>

namespace gmml
{
    namespace assembly
    {
        namespace
        {
            void insertNonIgnored(
                std::vector<size_t>& result,
                const std::vector<size_t>& indices,
                const std::vector<bool>& includeAtom,
                const std::vector<bool>& localIgnore)
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

            void insertIntersection(
                std::vector<size_t>& result,
                double overlapTolerance,
                const Sphere& sphere,
                const std::vector<Sphere>& bounds,
                const std::vector<bool>& includeAtom,
                const std::vector<size_t>& indices)
            {
                result.reserve(indices.size());
                for (size_t index : indices)
                {
                    if (includeAtom[index] && spheresOverlap(overlapTolerance, sphere, bounds[index]))
                    {
                        result.push_back(index);
                    }
                }
            }

            std::array<std::vector<size_t>, 2> atomsToCheck(
                const Graph& graph,
                const Bounds& bounds,
                const Selection& selectionA,
                const Selection& selectionB,
                const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
                double overlapTolerance,
                size_t residueA,
                size_t residueB)
            {
                const Sphere& residueBoundsA = bounds.residues[residueA];
                const Sphere& residueBoundsB = bounds.residues[residueB];
                const std::vector<size_t>& atomsA = residueAtoms(graph, residueA);
                const std::vector<size_t>& atomsB = residueAtoms(graph, residueB);
                std::vector<size_t> indicesA;
                std::vector<size_t> indicesB;
                size_t adjacency = util::indexOf(graph.residues.nodes.nodeAdjacencies[residueA], residueB);
                const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residueA];
                if (adjacency < edges.size())
                {
                    size_t edgeIndex = edges[adjacency];
                    const std::array<std::vector<bool>, 2>& ignoredAtoms = residueAtomsCloseToEdge[edgeIndex];
                    bool order = !(graph.residues.edges.nodeAdjacencies[edgeIndex][0] == residueA);
                    insertNonIgnored(indicesA, atomsA, selectionA.atoms, ignoredAtoms[order]);
                    insertNonIgnored(indicesB, atomsB, selectionB.atoms, ignoredAtoms[!order]);
                }
                else if (spheresOverlap(overlapTolerance, residueBoundsA, residueBoundsB))
                {
                    insertIntersection(
                        indicesA, overlapTolerance, residueBoundsB, bounds.atoms, selectionA.atoms, atomsA);
                    insertIntersection(
                        indicesB, overlapTolerance, residueBoundsA, bounds.atoms, selectionB.atoms, atomsB);
                }
                return {indicesA, indicesB};
            }
        } // namespace

        void insertIndicesOfIntersection(
            std::vector<size_t>& result,
            double overlapTolerance,
            const Sphere& sphere,
            const std::vector<Sphere>& coords,
            const std::vector<size_t>& indices)
        {
            result.reserve(indices.size());
            for (size_t index : indices)
            {
                auto& a = coords[index];
                if (spheresOverlap(overlapTolerance, sphere, a))
                {
                    result.push_back(index);
                }
            }
        }

        std::vector<size_t> intersectingIndices(
            double overlapTolerance,
            const Sphere& sphere,
            const std::vector<Sphere>& coords,
            const std::vector<size_t>& indices)
        {
            std::vector<size_t> result;
            insertIndicesOfIntersection(result, overlapTolerance, sphere, coords, indices);
            return result;
        }

        void addResidueOverlaps(
            std::vector<double>& result,
            const PotentialTable& potential,
            double overlapTolerance,
            const Graph& graph,
            const Bounds& bounds,
            const Selection& selectionA,
            const Selection& selectionB,
            const std::vector<Element>& atomElements,
            const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
            size_t residueA,
            size_t residueB)
        {
            std::array<std::vector<size_t>, 2> toCheck = atomsToCheck(
                graph, bounds, selectionA, selectionB, residueAtomsCloseToEdge, overlapTolerance, residueA, residueB);
            for (size_t atomA : toCheck[0])
            {
                Element elementA = atomElements[atomA];
                for (size_t atomB : toCheck[1])
                {
                    Element elementB = atomElements[atomB];
                    PotentialFactor factor = potentialFactor(potential, elementA, elementB);
                    double overlap = overlapAmount(factor, overlapTolerance, bounds.atoms[atomA], bounds.atoms[atomB]);
                    result[atomA] += overlap;
                    result[atomB] += overlap;
                }
            }
        }

        std::vector<double> overlapsBetweenSelections(
            const PotentialTable& potential,
            double overlapTolerance,
            const Graph& graph,
            const Bounds& bounds,
            const Selection& selectionA,
            const Selection& selectionB,
            const std::vector<Element>& atomElements,
            const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge)
        {
            std::vector<double> result(graph.indices.atomCount, 0.0);
            std::vector<size_t> moleculesA = selectedMolecules(selectionA);
            std::vector<size_t> moleculesB = selectedMolecules(selectionB);
            std::vector<std::vector<size_t>> moleculeResiduesA =
                moleculeSelectedResidues(graph, selectionA, moleculesA);
            std::vector<std::vector<size_t>> moleculeResiduesB =
                moleculeSelectedResidues(graph, selectionB, moleculesB);
            for (size_t maIndex = 0; maIndex < moleculesA.size(); maIndex++)
            {
                size_t moleculeA = moleculesA[maIndex];
                Sphere boundsA = bounds.molecules[moleculeA];
                for (size_t mbIndex = 0; mbIndex < moleculesB.size(); mbIndex++)
                {
                    size_t moleculeB = moleculesB[mbIndex];
                    Sphere boundsB = bounds.molecules[moleculeB];
                    if (spheresOverlap(overlapTolerance, bounds.molecules[moleculeA], bounds.molecules[moleculeB]))
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
                                    spheresOverlap(
                                        overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                                {
                                    addResidueOverlaps(
                                        result,
                                        potential,
                                        overlapTolerance,
                                        graph,
                                        bounds,
                                        selectionA,
                                        selectionB,
                                        atomElements,
                                        residueAtomsCloseToEdge,
                                        residueA,
                                        residueB);
                                }
                            }
                        }
                    }
                }
            }
            return result;
        }

        std::vector<double> overlapsWithinSelection(
            const PotentialTable& potential,
            double overlapTolerance,
            const Graph& graph,
            const Bounds& bounds,
            const Selection& selection,
            const std::vector<Element>& atomElements,
            const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge)
        {
            std::vector<double> result(graph.indices.atomCount, 0.0);
            std::vector<size_t> molecules = selectedMolecules(selection);
            std::vector<std::vector<size_t>> moleculeResidues = moleculeSelectedResidues(graph, selection, molecules);
            for (size_t maIndex = 0; maIndex < molecules.size(); maIndex++)
            {
                size_t moleculeA = molecules[maIndex];
                Sphere boundsA = bounds.molecules[moleculeA];
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
                                if (spheresOverlap(
                                        overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                                {
                                    addResidueOverlaps(
                                        result,
                                        potential,
                                        overlapTolerance,
                                        graph,
                                        bounds,
                                        selection,
                                        selection,
                                        atomElements,
                                        residueAtomsCloseToEdge,
                                        residueA,
                                        residueB);
                                }
                            }
                        }
                    }
                    else
                    {
                        Sphere boundsB = bounds.molecules[moleculeB];
                        std::vector<size_t> residuesA =
                            intersectingIndices(overlapTolerance, boundsB, bounds.residues, moleculeResidues[maIndex]);
                        std::vector<size_t> residuesB =
                            intersectingIndices(overlapTolerance, boundsA, bounds.residues, moleculeResidues[mbIndex]);
                        for (size_t residueA : residuesA)
                        {
                            for (size_t residueB : residuesB)
                            {
                                if (spheresOverlap(
                                        overlapTolerance, bounds.residues[residueA], bounds.residues[residueB]))
                                {
                                    addResidueOverlaps(
                                        result,
                                        potential,
                                        overlapTolerance,
                                        graph,
                                        bounds,
                                        selection,
                                        selection,
                                        atomElements,
                                        residueAtomsCloseToEdge,
                                        residueA,
                                        residueB);
                                }
                            }
                        }
                    }
                }
            }
            return result;
        }
    } // namespace assembly
} // namespace gmml
