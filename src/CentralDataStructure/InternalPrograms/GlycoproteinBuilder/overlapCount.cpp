#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    namespace
    {
        assembly::Selection selectMolecule(const assembly::Graph& graph, const assembly::Selection& selection, size_t n)
        {
            std::vector<bool> molecule(graph.indices.moleculeCount, false);
            molecule[n] = true;
            return assembly::selectByAtomsAndMolecules(graph, selection.atoms, molecule);
        };

        assembly::Selection selectGlycan(const assembly::Graph& graph, const AssemblyData& data,
                                         const assembly::Selection& selection, size_t n)
        {
            assembly::Selection glycanSelection = selectMolecule(graph, selection, data.glycans.moleculeId[n]);
            size_t residue                      = data.glycans.attachmentResidue[n];
            glycanSelection.residues[n]         = true;
            for (size_t atom : residueAtoms(graph, residue))
            {
                glycanSelection.atoms[atom] = selection.atoms[atom];
            }
            return glycanSelection;
        };
    } // namespace

    std::vector<cds::Overlap> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                            const assembly::Selection& selection, const assembly::Bounds& bounds)
    {
        return cds::overlapsWithinSelection(data.potentialTable, data.overlapTolerance, graph, bounds, selection,
                                            data.atoms.elements, data.residueEdges.atomsCloseToEdge);
    }

    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data,
                              const assembly::Selection& selection, const assembly::Bounds& bounds, size_t glycanId)
    {
        assembly::Selection glycanSelection               = selectGlycan(graph, data, selection, glycanId);
        std::vector<bool> otherMolecules                  = selection.molecules;
        otherMolecules[data.glycans.moleculeId[glycanId]] = false;
        assembly::Selection otherSelection =
            assembly::selectByAtomsAndMolecules(graph, selection.atoms, otherMolecules);

        std::vector<cds::Overlap> overlapWithin =
            cds::overlapsWithinSelection(data.potentialTable, data.overlapTolerance, graph, bounds, glycanSelection,
                                         data.atoms.elements, data.residueEdges.atomsCloseToEdge);
        std::vector<cds::Overlap> overlapBetween =
            cds::overlapsBetweenSelections(data.potentialTable, data.overlapTolerance, graph, bounds, glycanSelection,
                                           otherSelection, data.atoms.elements, data.residueEdges.atomsCloseToEdge);
        return cds::overlapVectorSum(overlapWithin) + cds::overlapVectorSum(overlapBetween);
    };

    std::vector<size_t> determineSitesWithOverlap(double threshold, const std::vector<size_t>& movedSites,
                                                  const assembly::Graph& graph, const AssemblyData& data,
                                                  const assembly::Selection& selection, const assembly::Bounds& bounds)

    {
        const std::vector<bool>& included = glycanIncluded(data, selection.molecules);

        auto hasProteinOverlap = [&](size_t n)
        {
            assembly::Selection glycanSelection  = selectGlycan(graph, data, selection, n);
            assembly::Selection proteinSelection = assembly::selectByAtomsAndMolecules(
                graph, selection.atoms,
                codeUtils::indicesToBools(graph.indices.moleculeCount, data.indices.proteinMolecules));

            return containsOverlapExceedingThreshold(
                threshold, cds::overlapsBetweenSelections(data.potentialTable, data.overlapTolerance, graph, bounds,
                                                          glycanSelection, proteinSelection, data.atoms.elements,
                                                          data.residueEdges.atomsCloseToEdge));
        };
        auto hasSelfOverlap = [&](size_t n)
        {
            std::vector<bool> glycanMolecule(graph.indices.moleculeCount, false);
            glycanMolecule[data.glycans.moleculeId[n]] = true;
            assembly::Selection glycanSelection =
                assembly::selectByAtomsAndMolecules(graph, selection.atoms, glycanMolecule);

            return containsOverlapExceedingThreshold(
                threshold,
                cds::overlapsWithinSelection(data.potentialTable, data.overlapTolerance, graph, bounds, glycanSelection,
                                             data.atoms.elements, data.residueEdges.atomsCloseToEdge));
        };
        size_t glycanCount = data.glycans.moleculeId.size();
        std::vector<size_t> sitesToCheck;
        sitesToCheck.reserve(glycanCount);
        for (size_t n : movedSites)
        {
            if (included[n])
            {
                sitesToCheck.push_back(n);
            }
        }
        std::vector<bool> justMoved(glycanCount, false);
        for (size_t n : sitesToCheck)
        {
            justMoved[n] = included[n];
        }
        std::vector<bool> glycanOverlap(glycanCount, false);
        for (size_t n : sitesToCheck)
        {
            for (size_t k = 0; k < glycanCount; k++)
            {
                bool avoidDoubleCount = k > n || !justMoved[k];
                if (included[k] && (k != n) && avoidDoubleCount && !(glycanOverlap[n] && glycanOverlap[k]))
                {
                    assembly::Selection selectionN    = selectGlycan(graph, data, selection, n);
                    assembly::Selection selectionK    = selectGlycan(graph, data, selection, k);
                    std::vector<cds::Overlap> overlap = cds::overlapsBetweenSelections(
                        data.potentialTable, data.overlapTolerance, graph, bounds, selectionN, selectionK,
                        data.atoms.elements, data.residueEdges.atomsCloseToEdge);
                    if (cds::containsOverlapExceedingThreshold(threshold, overlap))
                    {
                        glycanOverlap[n] = true;
                        glycanOverlap[k] = true;
                    }
                }
            }
        }
        std::vector<size_t> indices;
        for (size_t n = 0; n < glycanCount; n++)
        {
            // glycans which haven't moved won't overlap with protein or themselves (at least not more than before)
            if (included[n] && (glycanOverlap[n] || (justMoved[n] && (hasProteinOverlap(n) || hasSelfOverlap(n)))))
            {
                indices.push_back(n);
            }
        }
        return indices;
    }

} // namespace glycoproteinBuilder
