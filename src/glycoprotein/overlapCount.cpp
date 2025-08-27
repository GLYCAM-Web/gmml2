#include "include/glycoprotein/overlapCount.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyOverlap.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/glycoproteinUtil.hpp"
#include "include/util/containers.hpp"

#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            assembly::Selection selectMolecule(
                const assembly::Graph& graph, const assembly::Selection& selection, size_t n)
            {
                std::vector<bool> molecule(moleculeCount(graph.source), false);
                molecule[n] = true;
                return assembly::selectByAtomsAndMolecules(graph, selection.atoms, molecule);
            };

            assembly::Selection selectGlycan(
                const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection, size_t n)
            {
                assembly::Selection glycanSelection = selectMolecule(graph, selection, data.glycans.moleculeId[n]);
                size_t residue = data.glycans.attachmentResidue[n];
                glycanSelection.residues[n] = true;
                for (size_t atom : residueAtoms(graph, residue))
                {
                    glycanSelection.atoms[atom] = selection.atoms[atom];
                }
                return glycanSelection;
            };
        } // namespace

        std::vector<double> totalOverlaps(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds)
        {
            return overlapsWithinSelection(
                overlapSettings.potentialTable,
                overlapSettings.tolerance,
                graph,
                bounds,
                selection,
                data.atoms.elements,
                data.residueEdges.atomsCloseToEdge);
        }

        double localOverlap(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds,
            size_t glycanId)
        {
            assembly::Selection glycanSelection = selectGlycan(graph, data, selection, glycanId);
            std::vector<bool> otherMolecules = selection.molecules;
            otherMolecules[data.glycans.moleculeId[glycanId]] = false;
            assembly::Selection otherSelection =
                assembly::selectByAtomsAndMolecules(graph, selection.atoms, otherMolecules);

            std::vector<double> overlapWithin = overlapsWithinSelection(
                overlapSettings.potentialTable,
                overlapSettings.tolerance,
                graph,
                bounds,
                glycanSelection,
                data.atoms.elements,
                data.residueEdges.atomsCloseToEdge);
            std::vector<double> overlapBetween = overlapsBetweenSelections(
                overlapSettings.potentialTable,
                overlapSettings.tolerance,
                graph,
                bounds,
                glycanSelection,
                otherSelection,
                data.atoms.elements,
                data.residueEdges.atomsCloseToEdge);
            return overlapVectorSum(overlapWithin) + overlapVectorSum(overlapBetween);
        };

        OverlapSites determineOverlapState(
            double threshold,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const assembly::Selection& selection,
            const assembly::Bounds& bounds)
        {
            const std::vector<bool>& included = glycanIncluded(data, selection.molecules);

            size_t glycanCount = data.glycans.moleculeId.size();
            assembly::Selection proteinSelection =
                assembly::selectByAtomsAndMolecules(graph, selection.atoms, data.molecules.isProtein);
            std::vector<assembly::Selection> glycanSelections;
            glycanSelections.reserve(glycanCount);
            for (size_t n = 0; n < glycanCount; n++)
            {
                glycanSelections.push_back(selectGlycan(graph, data, selection, n));
            }

            auto hasProteinOverlap = [&](size_t n)
            {
                return containsOverlapExceedingThreshold(
                    threshold,
                    overlapsBetweenSelections(
                        overlapSettings.potentialTable,
                        overlapSettings.tolerance,
                        graph,
                        bounds,
                        glycanSelections[n],
                        proteinSelection,
                        data.atoms.elements,
                        data.residueEdges.atomsCloseToEdge));
            };
            auto hasSelfOverlap = [&](size_t n)
            {
                std::vector<bool> glycanMolecule(moleculeCount(graph.source), false);
                glycanMolecule[data.glycans.moleculeId[n]] = true;
                assembly::Selection glycanSelection =
                    assembly::selectByAtomsAndMolecules(graph, selection.atoms, glycanMolecule);

                return containsOverlapExceedingThreshold(
                    threshold,
                    overlapsWithinSelection(
                        overlapSettings.potentialTable,
                        overlapSettings.tolerance,
                        graph,
                        bounds,
                        glycanSelection,
                        data.atoms.elements,
                        data.residueEdges.atomsCloseToEdge));
            };
            std::vector<bool> hasGlycanOverlap(glycanCount, false);
            std::vector<std::vector<bool>> interactions(glycanCount, std::vector<bool>(glycanCount, false));
            for (size_t n = 0; n < glycanCount; n++)
            {
                for (size_t k = n + 1; k < glycanCount; k++)
                {
                    if (included[n] && included[k])
                    {
                        std::vector<double> overlap = overlapsBetweenSelections(
                            overlapSettings.potentialTable,
                            overlapSettings.tolerance,
                            graph,
                            bounds,
                            glycanSelections[n],
                            glycanSelections[k],
                            data.atoms.elements,
                            data.residueEdges.atomsCloseToEdge);
                        if (containsOverlapExceedingThreshold(threshold, overlap))
                        {
                            interactions[n][k] = true;
                            interactions[k][n] = true;
                            hasGlycanOverlap[n] = true;
                            hasGlycanOverlap[k] = true;
                        }
                    }
                }
            }
            std::vector<bool> hasOverlap(glycanCount, false);
            for (size_t n = 0; n < glycanCount; n++)
            {
                // glycans which haven't moved won't overlap with protein or themselves (at least not more than before)
                if (included[n] && (hasGlycanOverlap[n] || (hasProteinOverlap(n) || hasSelfOverlap(n))))
                {
                    hasOverlap[n] = true;
                }
            }

            std::vector<size_t> concertId = util::indexVector(glycanCount);
            for (size_t n = 0; n < glycanCount; n++)
            {
                for (size_t k = n + 1; k < glycanCount; k++)
                {
                    if (interactions[n][k])
                    {
                        concertId[k] = concertId[n];
                    }
                }
            }

            return {util::boolsToIndices(hasOverlap), hasOverlap, interactions, concertId};
        }
    } // namespace gpbuilder
} // namespace gmml
