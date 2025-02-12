#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    namespace
    {
        std::vector<bool> ignoredAtomsOf(const assembly::Graph& graph, size_t residueIndex, size_t atomIndex)
        {
            const std::vector<size_t>& atomAdj     = graph.atoms.nodes.nodeAdjacencies[atomIndex];
            const std::vector<size_t>& atomIndices = residueAtoms(graph, residueIndex);
            std::vector<bool> ignored(atomIndices.size(), false);
            ignored[codeUtils::indexOf(atomIndices, atomIndex)] = true;
            for (size_t adj : atomAdj)
            {
                if (graph.atomResidue[adj] == residueIndex)
                {
                    ignored[codeUtils::indexOf(atomIndices, adj)] = true;
                }
            }
            return ignored;
        };

        cds::BondedResidueOverlapInput bondedResidueOverlapInput(const assembly::Graph& graph, size_t bondIndex)
        {
            auto& atomBond             = graph.atoms.edges.nodeAdjacencies[bondIndex];
            size_t atomA               = atomBond[0];
            size_t atomB               = atomBond[1];
            size_t residueA            = graph.atomResidue[atomA];
            size_t residueB            = graph.atomResidue[atomB];
            std::vector<bool> ignoredA = ignoredAtomsOf(graph, residueA, atomA);
            std::vector<bool> ignoredB = ignoredAtomsOf(graph, residueB, atomB);
            return {
                {residueA, residueB},
                {ignoredA, ignoredB}
            };
        }

        std::vector<cds::BondedResidueOverlapInput> residueBonds(const assembly::Graph& graph, size_t residueA,
                                                                 size_t residueB)
        {
            std::vector<cds::BondedResidueOverlapInput> bonds;
            auto& adjacencies     = graph.residues.nodes.nodeAdjacencies[residueA];
            size_t adjacencyIndex = codeUtils::indexOf(adjacencies, residueB);
            if (adjacencyIndex < adjacencies.size())
            {
                size_t edgeIndex = graph.residues.nodes.edgeAdjacencies[residueA][adjacencyIndex];
                size_t bondIndex = residueEdgeToAtomEdgeIndex(graph, edgeIndex);
                bonds.push_back(bondedResidueOverlapInput(graph, bondIndex));
            }
            return bonds;
        }

        std::vector<cds::BondedResidueOverlapInput> moleculeBonds(const assembly::Graph& graph, size_t moleculeA,
                                                                  size_t moleculeB)
        {
            std::vector<cds::BondedResidueOverlapInput> bonds;
            auto& adjacencies     = graph.molecules.nodes.nodeAdjacencies[moleculeA];
            size_t adjacencyIndex = codeUtils::indexOf(adjacencies, moleculeB);
            if (adjacencyIndex < adjacencies.size())
            {
                size_t edgeIndex     = graph.molecules.nodes.edgeAdjacencies[moleculeA][adjacencyIndex];
                size_t atomBondIndex = moleculeEdgeToAtomEdgeIndex(graph, edgeIndex);
                bonds.push_back(bondedResidueOverlapInput(graph, atomBondIndex));
            }
            return bonds;
        }
    } // namespace

    std::vector<cds::Overlap> intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                  const MutableData& mutableData,
                                                  const std::vector<double>& residueOverlapWeights,
                                                  const std::vector<bool>& includedAtoms, size_t glycanId)
    {
        const std::vector<size_t>& glycanLinkages = data.glycans.linkages[glycanId];
        std::vector<cds::Overlap> result(graph.residueCount, {0.0, 0.0});
        // skip first linkage as it connects to protein. We're only counting glycan atoms here
        for (size_t n = 1; n < glycanLinkages.size(); n++)
        {
            const ResidueLinkageIndices& linkage = data.indices.residueLinkages[glycanLinkages[n]];
            // only take first non-reducing residue to avoid any double-counting
            size_t residueA                      = linkage.nonReducingResidues[0];
            const std::vector<size_t>& residuesB = linkage.reducingResidues;

            std::vector<cds::BondedResidueOverlapInput> bonds = residueBonds(graph, residueA, residuesB[0]);
            cds::addOverlapsTo(result, cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                                                  graph.residues.nodes.elements, residueOverlapWeights,
                                                                  includedAtoms, bonds, {residueA}, residuesB));
        }
        return result;
    }

    std::vector<cds::Overlap> moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                               const MutableData& mutableData,
                                               const std::vector<double>& residueOverlapWeights,
                                               const std::vector<bool>& includedAtoms, size_t moleculeA,
                                               size_t moleculeB)
    {
        const cds::Sphere& boundsA = mutableData.moleculeBounds[moleculeA];
        const cds::Sphere& boundsB = mutableData.moleculeBounds[moleculeB];
        if (!cds::spheresOverlap(constants::overlapTolerance, boundsA, boundsB))
        {
            return std::vector<cds::Overlap>(graph.residueCount, {0.0, 0.0});
        }
        else
        {
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graph, moleculeA, moleculeB);
            std::vector<size_t> residuesA =
                cds::intersectingIndices(boundsB, mutableData.residueBounds, moleculeResidues(graph, moleculeA));
            std::vector<size_t> residuesB =
                cds::intersectingIndices(boundsA, mutableData.residueBounds, moleculeResidues(graph, moleculeB));
            return cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                              graph.residues.nodes.elements, residueOverlapWeights, includedAtoms,
                                              bonds, residuesA, residuesB);
        }
    }

    std::vector<cds::Overlap> moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                                      const MutableData& mutableData,
                                                      const std::vector<bool>& includedAtoms, size_t molecule,
                                                      size_t residue)
    {
        const cds::Sphere& moleculeBounds = mutableData.moleculeBounds[molecule];
        const cds::Sphere& residueBounds  = mutableData.residueBounds[residue];
        if (!cds::spheresOverlap(constants::overlapTolerance, moleculeBounds, residueBounds))
        {
            return std::vector<cds::Overlap>(graph.residueCount, {0.0, 0.0});
        }
        else
        {
            size_t residueMolecule                            = graph.residueMolecule[residue];
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graph, molecule, residueMolecule);
            return cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                              graph.residues.nodes.elements, data.residues.overlapWeights,
                                              includedAtoms, bonds, {residue}, moleculeResidues(graph, molecule));
        }
    }

    std::vector<cds::Overlap> totalOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                            const MutableData& mutableData,
                                            const std::vector<double>& residueOverlapWeights,
                                            const std::vector<bool>& includedAtoms, OverlapWeight weight)
    {
        std::vector<cds::Overlap> result(graph.residueCount, {0.0, 0.0});
        const std::vector<size_t> glycans = includedGlycanIndices(data, mutableData);

        for (size_t n : glycans)
        {
            size_t glycanMolecule = data.glycans.moleculeId[n];
            cds::addOverlapsTo(
                result, cds::scaledOverlaps(weight.self, intraGlycanOverlaps(graph, data, mutableData,
                                                                             residueOverlapWeights, includedAtoms, n)));
            for (size_t k : data.indices.proteinMolecules)
            {
                cds::addOverlapsTo(result, moleculeOverlaps(graph, data, mutableData, residueOverlapWeights,
                                                            includedAtoms, glycanMolecule, k));
            }
            for (size_t k : glycans)
            {
                if (k > n)
                {
                    cds::addOverlapsTo(result,
                                       moleculeOverlaps(graph, data, mutableData, residueOverlapWeights, includedAtoms,
                                                        glycanMolecule, data.glycans.moleculeId[k]));
                }
            }
        }

        return result;
    }

    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                              const std::vector<bool>& includedAtoms, size_t glycanId, double selfWeight)
    {
        const std::vector<double>& residueOverlapWeights = data.residues.overlapWeights;
        const std::vector<size_t> glycans                = includedGlycanIndices(data, mutableData);
        cds::Overlap overlap                             = cds::overlapVectorSum(intraGlycanOverlaps(
                                   graph, data, mutableData, residueOverlapWeights, includedAtoms, glycanId)) *
                               selfWeight;
        size_t glycanMolecule = data.glycans.moleculeId[glycanId];
        for (size_t n : data.indices.proteinMolecules)
        {
            overlap += cds::overlapVectorSum(
                moleculeOverlaps(graph, data, mutableData, residueOverlapWeights, includedAtoms, n, glycanMolecule));
        }

        for (size_t n : glycans)
        {
            if (n != glycanId)
            {
                size_t other = data.glycans.moleculeId[n];
                overlap      += cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData, residueOverlapWeights,
                                                                       includedAtoms, glycanMolecule, other));
                overlap += cds::overlapVectorSum(moleculeResidueOverlaps(graph, data, mutableData, includedAtoms, other,
                                                                         data.glycans.attachmentResidue[glycanId]));
            }
        }
        return overlap;
    };

    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const MutableData& mutableData,
                                                  const std::vector<bool>& includedAtoms)
    {
        const std::vector<double>& residueOverlapWeights = data.residues.overlapWeights;
        const std::vector<bool>& included                = glycanIncluded(data, mutableData);
        auto hasProteinOverlap                           = [&](size_t n)
        {
            cds::Overlap overlap {0.0, 0.0};
            for (size_t k : data.indices.proteinMolecules)
            {
                overlap += cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData, residueOverlapWeights,
                                                                  includedAtoms, k, data.glycans.moleculeId[n]));
            }
            return overlap.count > 0;
        };
        auto hasSelfOverlap = [&](size_t n)
        {
            cds::Overlap overlap = cds::overlapVectorSum(
                intraGlycanOverlaps(graph, data, mutableData, residueOverlapWeights, includedAtoms, n));
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
                    if (cds::overlapVectorSum(moleculeOverlaps(graph, data, mutableData, residueOverlapWeights,
                                                               includedAtoms, data.glycans.moleculeId[n],
                                                               data.glycans.moleculeId[k]))
                                .count > 0 ||
                        cds::overlapVectorSum(moleculeResidueOverlaps(graph, data, mutableData, includedAtoms,
                                                                      data.glycans.moleculeId[k],
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
