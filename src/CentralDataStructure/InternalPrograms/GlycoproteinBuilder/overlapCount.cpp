#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::Overlap> intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                  const assembly::Bounds& bounds,
                                                  const std::vector<double>& residueWeights,
                                                  const std::vector<bool>& includedAtoms, size_t glycanId)
    {
        const std::vector<size_t>& glycanLinkages = data.glycans.linkages[glycanId];
        std::vector<cds::Overlap> result(graph.atomCount, {0.0, 0.0});
        // skip first linkage as it connects to protein. We're only counting glycan atoms here
        for (size_t n = 1; n < glycanLinkages.size(); n++)
        {
            const ResidueLinkageIndices& linkage = data.indices.residueLinkages[glycanLinkages[n]];
            // only take first non-reducing residue to avoid any double-counting
            size_t residueA                      = linkage.nonReducingResidues[0];
            const std::vector<size_t>& residuesB = linkage.reducingResidues;
            cds::addOverlapsTo(result,
                               cds::CountOverlappingAtoms(data.potentialTable, data.overlapProperties, graph,
                                                          {bounds.atoms, data.atoms.elements, includedAtoms},
                                                          {bounds.residues, residueWeights},
                                                          data.residueEdges.atomsCloseToEdge, {residueA}, residuesB));
        }
        return result;
    }

    std::vector<cds::Overlap> moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                               const assembly::Bounds& bounds,
                                               const std::vector<double>& residueWeights,
                                               const std::vector<bool>& includedAtoms, size_t moleculeA,
                                               size_t moleculeB)
    {
        double overlapTolerance    = data.overlapProperties.tolerance;
        const cds::Sphere& boundsA = bounds.molecules[moleculeA];
        const cds::Sphere& boundsB = bounds.molecules[moleculeB];
        if (!cds::spheresOverlap(overlapTolerance, boundsA, boundsB))
        {
            return std::vector<cds::Overlap>(graph.atomCount, {0.0, 0.0});
        }
        else
        {
            std::vector<size_t> residuesA = cds::intersectingIndices(overlapTolerance, boundsB, bounds.residues,
                                                                     moleculeResidues(graph, moleculeA));
            std::vector<size_t> residuesB = cds::intersectingIndices(overlapTolerance, boundsA, bounds.residues,
                                                                     moleculeResidues(graph, moleculeB));
            return cds::CountOverlappingAtoms(
                data.potentialTable, data.overlapProperties, graph, {bounds.atoms, data.atoms.elements, includedAtoms},
                {bounds.residues, residueWeights}, data.residueEdges.atomsCloseToEdge, residuesA, residuesB);
        }
    }

    std::vector<cds::Overlap> moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                      const assembly::Bounds& bounds,
                                                      const std::vector<double>& residueWeights,
                                                      const std::vector<bool>& includedAtoms, size_t molecule,
                                                      size_t residue)
    {
        const cds::Sphere& moleculeBounds = bounds.molecules[molecule];
        const cds::Sphere& residueBounds  = bounds.residues[residue];
        if (!cds::spheresOverlap(data.overlapProperties.tolerance, moleculeBounds, residueBounds))
        {
            return std::vector<cds::Overlap>(graph.atomCount, {0.0, 0.0});
        }
        else
        {
            return cds::CountOverlappingAtoms(data.potentialTable, data.overlapProperties, graph,
                                              {bounds.atoms, data.atoms.elements, includedAtoms},
                                              {bounds.residues, residueWeights}, data.residueEdges.atomsCloseToEdge,
                                              {residue}, moleculeResidues(graph, molecule));
        }
    }

    std::vector<cds::Overlap> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                            const MutableData& mutableData, const std::vector<double>& residueWeights,
                                            const std::vector<bool>& includedAtoms, OverlapMultiplier overlapMultiplier)
    {
        std::vector<cds::Overlap> result(graph.atomCount, {0.0, 0.0});
        const std::vector<size_t> glycans = includedGlycanIndices(data, mutableData);

        for (size_t n : glycans)
        {
            size_t glycanMolecule = data.glycans.moleculeId[n];
            cds::addOverlapsTo(result, cds::scaledOverlaps(overlapMultiplier.self,
                                                           intraGlycanOverlaps(graph, data, mutableData.bounds,
                                                                               residueWeights, includedAtoms, n)));
            for (size_t k : data.indices.proteinMolecules)
            {
                cds::addOverlapsTo(result, moleculeOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                            includedAtoms, glycanMolecule, k));
            }
            for (size_t k : glycans)
            {
                if (k > n)
                {
                    cds::addOverlapsTo(result,
                                       moleculeOverlaps(graph, data, mutableData.bounds, residueWeights, includedAtoms,
                                                        glycanMolecule, data.glycans.moleculeId[k]));
                }
            }
        }

        return result;
    }

    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                              const std::vector<double>& residueWeights, const std::vector<bool>& includedAtoms,
                              size_t glycanId, double selfWeight)
    {
        const std::vector<size_t> glycans = includedGlycanIndices(data, mutableData);
        cds::Overlap overlap              = cds::overlapVectorSum(intraGlycanOverlaps(graph, data, mutableData.bounds,
                                                                                      residueWeights, includedAtoms, glycanId)) *
                               selfWeight;
        size_t glycanMolecule = data.glycans.moleculeId[glycanId];
        for (size_t n : data.indices.proteinMolecules)
        {
            overlap += cds::overlapVectorSum(
                moleculeOverlaps(graph, data, mutableData.bounds, residueWeights, includedAtoms, n, glycanMolecule));
        }

        for (size_t n : glycans)
        {
            if (n != glycanId)
            {
                size_t other = data.glycans.moleculeId[n];
                overlap      += cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                                       includedAtoms, glycanMolecule, other));
                overlap      += cds::overlapVectorSum(moleculeResidueOverlaps(graph, data, mutableData.bounds,
                                                                              residueWeights, includedAtoms, other,
                                                                              data.glycans.attachmentResidue[glycanId]));
            }
        }
        return overlap;
    };

    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const MutableData& mutableData,
                                                  const std::vector<bool>& includedAtoms)
    {
        const std::vector<double>& residueWeights = data.equalResidueWeight;
        const std::vector<bool>& included         = glycanIncluded(data, mutableData);
        auto hasProteinOverlap                    = [&](size_t n)
        {
            cds::Overlap overlap {0.0, 0.0};
            for (size_t k : data.indices.proteinMolecules)
            {
                overlap += cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                                  includedAtoms, k, data.glycans.moleculeId[n]));
            }
            return overlap.count > 0;
        };
        auto hasSelfOverlap = [&](size_t n)
        {
            cds::Overlap overlap = cds::overlapVectorSum(
                intraGlycanOverlaps(graph, data, mutableData.bounds, residueWeights, includedAtoms, n));
            return overlap.count > 0;
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
                    if (cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                               includedAtoms, data.glycans.moleculeId[n],
                                                               data.glycans.moleculeId[k]))
                                .count > 0 ||
                        cds::overlapVectorSum(moleculeResidueOverlaps(graph, data, mutableData.bounds, residueWeights,
                                                                      includedAtoms, data.glycans.moleculeId[k],
                                                                      data.glycans.attachmentResidue[n]))
                                .count > 0)
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
