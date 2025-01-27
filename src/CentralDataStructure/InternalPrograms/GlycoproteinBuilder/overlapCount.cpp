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

    cds::Overlap intraGlycanOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                     const MutableData& mutableData, size_t glycanId)
    {
        const std::vector<size_t>& glycanLinkages = data.indices.glycans[glycanId].linkages;
        cds::Overlap overlap                      = {0.0, 0.0};
        // skip first linkage as it connects to protein. We're only counting glycan atoms here
        for (size_t n = 1; n < glycanLinkages.size(); n++)
        {
            const ResidueLinkageIndices& linkage = data.indices.residueLinkages[glycanLinkages[n]];
            // only take first non-reducing residue to avoid any double-counting
            size_t residueA                      = linkage.nonReducingResidues[0];
            const std::vector<size_t>& residuesB = linkage.reducingResidues;

            std::vector<cds::BondedResidueOverlapInput> bonds = residueBonds(graph, residueA, residuesB[0]);
            overlap += cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                                  graph.residues.nodes.elements, data.residues.overlapWeights,
                                                  mutableData.atomIgnored, bonds, {residueA}, residuesB);
        }
        return overlap;
    }

    cds::Overlap moleculeOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                  const MutableData& mutableData, size_t moleculeA, size_t moleculeB)
    {
        const cds::Sphere& boundsA = mutableData.moleculeBounds[moleculeA];
        const cds::Sphere& boundsB = mutableData.moleculeBounds[moleculeB];
        if (!cds::spheresOverlap(constants::overlapTolerance, boundsA, boundsB))
        {
            return cds::Overlap {0.0, 0.0};
        }
        else
        {
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graph, moleculeA, moleculeB);
            std::vector<size_t> residuesA;
            std::vector<size_t> residuesB;
            cds::insertIndicesOfIntersection(residuesA, boundsB, mutableData.residueBounds,
                                             moleculeResidues(graph, moleculeA));
            cds::insertIndicesOfIntersection(residuesB, boundsA, mutableData.residueBounds,
                                             moleculeResidues(graph, moleculeB));
            return cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                              graph.residues.nodes.elements, data.residues.overlapWeights,
                                              mutableData.atomIgnored, bonds, residuesA, residuesB);
        }
    }

    cds::Overlap moleculeResidueOverlaps(const assembly::Graph& graph, const AssemblyData& data,
                                         const MutableData& mutableData, size_t molecule, size_t residue)
    {
        const cds::Sphere& moleculeBounds = mutableData.moleculeBounds[molecule];
        const cds::Sphere& residueBounds  = mutableData.residueBounds[residue];
        if (!cds::spheresOverlap(constants::overlapTolerance, moleculeBounds, residueBounds))
        {
            return cds::Overlap {0.0, 0.0};
        }
        else
        {
            size_t residueMolecule                            = graph.residueMolecule[residue];
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graph, molecule, residueMolecule);
            return cds::CountOverlappingAtoms(mutableData.atomBounds, mutableData.residueBounds,
                                              graph.residues.nodes.elements, data.residues.overlapWeights,
                                              mutableData.atomIgnored, bonds, {residue},
                                              moleculeResidues(graph, molecule));
        }
    }

    cds::Overlap totalOverlaps(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                               OverlapWeight weight)
    {
        cds::Overlap overlap {0.0, 0.0};
        const std::vector<GlycanIndices>& glycans = data.indices.glycans;
        const std::vector<bool>& included         = mutableData.glycanIncluded;

        for (size_t n = 0; n < glycans.size(); n++)
        {
            if (included[n])
            {
                size_t glycanMolecule = glycans[n].glycanMolecule;
                overlap               += intraGlycanOverlaps(graph, data, mutableData, n) * weight.self;
                for (size_t k : data.indices.proteinMolecules)
                {
                    overlap += moleculeOverlaps(graph, data, mutableData, glycanMolecule, k);
                }
                for (size_t k = n + 1; k < glycans.size(); k++)
                {
                    if (included[k])
                    {
                        overlap +=
                            moleculeOverlaps(graph, data, mutableData, glycanMolecule, glycans[k].glycanMolecule);
                    }
                }
            }
        }

        return overlap;
    }

    cds::Overlap localOverlap(const assembly::Graph& graph, const AssemblyData& data, const MutableData& mutableData,
                              size_t glycanId, double selfWeight)
    {
        const std::vector<GlycanIndices>& glycans = data.indices.glycans;
        const GlycanIndices& thisGlycan           = glycans[glycanId];
        cds::Overlap overlap = intraGlycanOverlaps(graph, data, mutableData, glycanId) * selfWeight;
        for (size_t n : data.indices.proteinMolecules)
        {
            overlap += moleculeOverlaps(graph, data, mutableData, n, thisGlycan.glycanMolecule);
        }

        for (size_t n = 0; n < glycans.size(); n++)
        {
            if (mutableData.glycanIncluded[n] && n != glycanId)
            {
                size_t other = glycans[n].glycanMolecule;
                overlap      += moleculeOverlaps(graph, data, mutableData, thisGlycan.glycanMolecule, other);
                overlap      += moleculeResidueOverlaps(graph, data, mutableData, other, thisGlycan.attachmentResidue);
            }
        }
        return overlap;
    };

    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const assembly::Graph& graph,
                                                  const AssemblyData& data, const MutableData& mutableData)
    {
        const std::vector<GlycanIndices>& glycans = data.indices.glycans;
        const std::vector<bool>& included         = mutableData.glycanIncluded;
        auto hasProteinOverlap                    = [&](size_t n)
        {
            cds::Overlap overlap {0.0, 0.0};
            for (size_t k : data.indices.proteinMolecules)
            {
                overlap += moleculeOverlaps(graph, data, mutableData, k, glycans[n].glycanMolecule);
            }
            return overlap.count > 0;
        };
        auto hasSelfOverlap = [&](size_t n)
        {
            auto overlap = intraGlycanOverlaps(graph, data, mutableData, n);
            return overlap.count > 0;
        };
        std::vector<size_t> sitesToCheck;
        sitesToCheck.reserve(glycans.size());
        for (size_t n : movedSites)
        {
            if (included[n])
            {
                sitesToCheck.push_back(n);
            }
        }
        std::vector<bool> justMoved(glycans.size(), false);
        for (size_t n : sitesToCheck)
        {
            justMoved[n] = included[n];
        }
        std::vector<bool> glycanOverlap(glycans.size(), false);
        for (size_t n : sitesToCheck)
        {
            for (size_t k = 0; k < glycans.size(); k++)
            {
                bool avoidDoubleCount = k > n || !justMoved[k];
                if (included[k] && (k != n) && avoidDoubleCount && !(glycanOverlap[n] && glycanOverlap[k]))
                {
                    if (moleculeOverlaps(graph, data, mutableData, glycans[n].glycanMolecule, glycans[k].glycanMolecule)
                                .count > 0 ||
                        moleculeResidueOverlaps(graph, data, mutableData, glycans[k].glycanMolecule,
                                                glycans[n].attachmentResidue)
                                .count > 0)
                    {
                        glycanOverlap[n] = true;
                        glycanOverlap[k] = true;
                    }
                }
            }
        }
        std::vector<size_t> indices;
        for (size_t n = 0; n < glycans.size(); n++)
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
