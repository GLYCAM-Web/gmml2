#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    namespace
    {
        std::vector<bool> ignoredAtomsOf(const AssemblyGraphs& graphs, size_t residueIndex, size_t atomIndex)
        {
            const std::vector<size_t>& atomAdj     = graphs.atoms.nodes.nodeAdjacencies[atomIndex];
            const std::vector<size_t>& atomIndices = residueAtoms(graphs, residueIndex);
            std::vector<bool> ignored(atomIndices.size(), false);
            ignored[codeUtils::indexOf(atomIndices, atomIndex)] = true;
            for (size_t adj : atomAdj)
            {
                if (graphs.indices.atomResidue[adj] == residueIndex)
                {
                    ignored[codeUtils::indexOf(atomIndices, adj)] = true;
                }
            }
            return ignored;
        };

        cds::BondedResidueOverlapInput bondedResidueOverlapInput(const AssemblyGraphs& graphs, size_t bondIndex)
        {
            auto& atomBond             = graphs.atoms.edges.nodeAdjacencies[bondIndex];
            size_t atomA               = atomBond[0];
            size_t atomB               = atomBond[1];
            size_t residueA            = graphs.indices.atomResidue[atomA];
            size_t residueB            = graphs.indices.atomResidue[atomB];
            std::vector<bool> ignoredA = ignoredAtomsOf(graphs, residueA, atomA);
            std::vector<bool> ignoredB = ignoredAtomsOf(graphs, residueB, atomB);
            return {
                {residueA, residueB},
                {ignoredA, ignoredB}
            };
        }

        std::vector<cds::BondedResidueOverlapInput> residueBonds(const AssemblyGraphs& graphs, size_t residueA,
                                                                 size_t residueB)
        {
            std::vector<cds::BondedResidueOverlapInput> bonds;
            auto& adjacencies     = graphs.residues.nodes.nodeAdjacencies[residueA];
            size_t adjacencyIndex = codeUtils::indexOf(adjacencies, residueB);
            if (adjacencyIndex < adjacencies.size())
            {
                size_t edgeIndex = graphs.residues.nodes.edgeAdjacencies[residueA][adjacencyIndex];
                size_t bondIndex = graphs.residues.edges.indices[edgeIndex];
                bonds.push_back(bondedResidueOverlapInput(graphs, bondIndex));
            }
            return bonds;
        }

        std::vector<cds::BondedResidueOverlapInput> moleculeBonds(const AssemblyGraphs& graphs, size_t moleculeA,
                                                                  size_t moleculeB)
        {
            std::vector<cds::BondedResidueOverlapInput> bonds;
            auto& adjacencies     = graphs.molecules.nodes.nodeAdjacencies[moleculeA];
            size_t adjacencyIndex = codeUtils::indexOf(adjacencies, moleculeB);
            if (adjacencyIndex < adjacencies.size())
            {
                size_t edgeIndex     = graphs.molecules.nodes.edgeAdjacencies[moleculeA][adjacencyIndex];
                size_t atomBondIndex = moleculeEdgeToAtomEdgeIndex(graphs, edgeIndex);
                bonds.push_back(bondedResidueOverlapInput(graphs, atomBondIndex));
            }
            return bonds;
        }
    } // namespace

    cds::Overlap intraGlycanOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId)
    {
        auto& glycanLinkages = graphs.glycans[glycanId].linkages;
        cds::Overlap overlap = {0.0, 0.0};
        // skip first linkage as it connects to protein. We're only counting glycan atoms here
        for (size_t n = 1; n < glycanLinkages.size(); n++)
        {
            auto& linkage                        = graphs.residueLinkages[glycanLinkages[n]];
            // only take first non-reducing residue to avoid any double-counting
            size_t residueA                      = linkage.nonReducingResidues[0];
            const std::vector<size_t>& residuesB = linkage.reducingResidues;

            std::vector<cds::BondedResidueOverlapInput> bonds = residueBonds(graphs, residueA, residuesB[0]);
            overlap +=
                cds::CountOverlappingAtoms(data.atoms.bounds, data.residues.bounds, graphs.residues.nodes.elements,
                                           data.residues.overlapWeights, bonds, {residueA}, residuesB);
        }
        return overlap;
    }

    cds::Overlap moleculeOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t moleculeA,
                                  size_t moleculeB)
    {
        const cds::Sphere& boundsA = data.molecules.bounds[moleculeA];
        const cds::Sphere& boundsB = data.molecules.bounds[moleculeB];
        if (!cds::spheresOverlap(constants::overlapTolerance, boundsA, boundsB))
        {
            return cds::Overlap {0.0, 0.0};
        }
        else
        {
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graphs, moleculeA, moleculeB);
            std::vector<size_t> residuesA;
            std::vector<size_t> residuesB;
            cds::insertIndicesOfIntersection(residuesA, boundsB, data.residues.bounds,
                                             moleculeResidues(graphs, moleculeA));
            cds::insertIndicesOfIntersection(residuesB, boundsA, data.residues.bounds,
                                             moleculeResidues(graphs, moleculeB));
            return cds::CountOverlappingAtoms(data.atoms.bounds, data.residues.bounds, graphs.residues.nodes.elements,
                                              data.residues.overlapWeights, bonds, residuesA, residuesB);
        }
    }

    cds::Overlap moleculeResidueOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t molecule,
                                         size_t residue)
    {
        const cds::Sphere& moleculeBounds = data.molecules.bounds[molecule];
        const cds::Sphere& residueBounds  = data.residues.bounds[residue];
        if (!cds::spheresOverlap(constants::overlapTolerance, moleculeBounds, residueBounds))
        {
            return cds::Overlap {0.0, 0.0};
        }
        else
        {
            size_t residueMolecule                            = graphs.indices.residueMolecule[residue];
            std::vector<cds::BondedResidueOverlapInput> bonds = moleculeBonds(graphs, molecule, residueMolecule);
            return cds::CountOverlappingAtoms(data.atoms.bounds, data.residues.bounds, graphs.residues.nodes.elements,
                                              data.residues.overlapWeights, bonds, {residue},
                                              moleculeResidues(graphs, molecule));
        }
    }

    cds::Overlap totalOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, OverlapWeight weight)
    {
        cds::Overlap overlap {0.0, 0.0};
        const std::vector<GlycanIndices>& glycans = graphs.glycans;
        const std::vector<bool>& included         = data.glycanData.included;

        for (size_t n = 0; n < glycans.size(); n++)
        {
            if (included[n])
            {
                size_t glycanMolecule = glycans[n].glycanMolecule;
                overlap               += intraGlycanOverlaps(graphs, data, n) * weight.self;
                for (size_t k : graphs.proteinMolecules)
                {
                    overlap += moleculeOverlaps(graphs, data, glycanMolecule, k);
                }
                for (size_t k = n + 1; k < glycans.size(); k++)
                {
                    if (included[k])
                    {
                        overlap += moleculeOverlaps(graphs, data, glycanMolecule, glycans[k].glycanMolecule);
                    }
                }
            }
        }

        return overlap;
    }

    cds::Overlap localOverlap(const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId,
                              double selfWeight)
    {
        const std::vector<GlycanIndices>& glycans = graphs.glycans;
        const GlycanIndices& thisGlycan           = glycans[glycanId];
        cds::Overlap overlap                      = intraGlycanOverlaps(graphs, data, glycanId) * selfWeight;
        for (size_t n : graphs.proteinMolecules)
        {
            overlap += moleculeOverlaps(graphs, data, n, thisGlycan.glycanMolecule);
        }

        for (size_t n = 0; n < glycans.size(); n++)
        {
            if (data.glycanData.included[n] && n != glycanId)
            {
                size_t other = glycans[n].glycanMolecule;
                overlap      += moleculeOverlaps(graphs, data, thisGlycan.glycanMolecule, other);
                overlap      += moleculeResidueOverlaps(graphs, data, other, thisGlycan.attachmentResidue);
            }
        }
        return overlap;
    };

    std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const AssemblyGraphs& graphs,
                                                  const AssemblyData& data)
    {
        const std::vector<GlycanIndices>& glycans = graphs.glycans;
        const std::vector<bool>& included         = data.glycanData.included;
        auto hasProteinOverlap                    = [&](size_t n)
        {
            cds::Overlap overlap {0.0, 0.0};
            for (size_t k : graphs.proteinMolecules)
            {
                overlap += moleculeOverlaps(graphs, data, k, glycans[n].glycanMolecule);
            }
            return overlap.count > 0;
        };
        auto hasSelfOverlap = [&](size_t n)
        {
            auto overlap = intraGlycanOverlaps(graphs, data, n);
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
                    if (moleculeOverlaps(graphs, data, glycans[n].glycanMolecule, glycans[k].glycanMolecule).count >
                            0 ||
                        moleculeResidueOverlaps(graphs, data, glycans[k].glycanMolecule, glycans[n].attachmentResidue)
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
