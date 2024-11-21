#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"

#include <array>
#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly toGlycoproteinAssemblyStructs(std::vector<cds::Molecule*>& molecules,
                                                       std::vector<GlycosylationSite>& glycosites,
                                                       const OverlapWeight overlapWeight)
    {
        cds::GraphIndexData graphIndices = cds::toIndexData(molecules);
        graph::Database atomGraphData    = cds::createGraphData(graphIndices);
        graph::Graph atomGraph           = graph::identity(atomGraphData);
        graph::Graph residueGraph        = graph::quotient(atomGraphData, graphIndices.atomResidue);
        graph::Graph moleculeGraph       = graph::quotient(graph::asData(residueGraph), graphIndices.residueMolecule);
        std::vector<cds::Atom*>& atoms   = graphIndices.atoms;
        std::vector<cds::Residue*>& residues = graphIndices.residues;

        std::vector<cds::Sphere> atomBoundingSpheres = cds::atomCoordinatesWithRadii(atoms);
        AtomData atomData {
            cds::atomNames(atoms),    cds::atomTypes(atoms),   cds::atomNumbers(atoms), cds::atomAtomicNumbers(atoms),
            cds::atomElements(atoms), cds::atomCharges(atoms), atomBoundingSpheres};

        auto boundingSpheresOf =
            [](const std::vector<cds::Sphere>& spheres, const std::vector<std::vector<size_t>>& indexVector)
        {
            std::vector<cds::Sphere> result;
            result.reserve(indexVector.size());
            for (auto& indices : indexVector)
            {
                result.push_back(cds::boundingSphere(codeUtils::indexValues(spheres, indices)));
            }
            return result;
        };

        auto indexOfResidues = [&residues](const std::vector<cds::Residue*>& toFind)
        {
            std::vector<size_t> result;
            result.reserve(toFind.size());
            for (auto& res : toFind)
            {
                result.push_back(codeUtils::indexOf(residues, res));
            }
            return result;
        };

        auto indexOfAtoms = [&atoms](const std::vector<cds::Atom*>& toFind)
        {
            std::vector<size_t> result;
            result.reserve(toFind.size());
            for (auto& res : toFind)
            {
                result.push_back(codeUtils::indexOf(atoms, res));
            }
            return result;
        };

        auto closelyBondedAtoms = [&](size_t residueId, size_t atomId)
        {
            if (graphIndices.atomResidue[atomId] != residueId)
            {
                throw std::runtime_error("bork");
            }
            const std::vector<size_t>& elements = residueGraph.nodes.elements[residueId];
            std::vector<bool> result(elements.size(), false);
            for (size_t n = 0; n < elements.size(); n++)
            {
                const std::vector<size_t>& adjacencies = atomGraph.nodes.nodeAdjacencies[atomId];
                size_t id                              = elements[n];
                result[n]                              = (id == atomId) || codeUtils::contains(adjacencies, id);
            }
            return result;
        };

        std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
        for (auto& a : glycosites)
        {
            glycosidicLinkages.push_back(a.GetGlycan()->GetGlycosidicLinkages());
        }

        std::vector<MoleculeType> moleculeTypes(graphIndices.molecules.size(), MoleculeType::protein);
        std::vector<cds::AngleWithMetadata> rotatableDihedralCurrentShape;
        std::vector<RotatableDihedralIndices> rotatableDihedralIndices;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycamMetadata::RotamerType> linkageRotamerTypes;
        std::vector<std::vector<GlycamMetadata::DihedralAngleDataVector>> linkageMetadata;
        std::vector<std::vector<cds::BondedResidueOverlapInput>> linkageOverlapBonds;
        std::vector<bool> linkageBranching;
        std::vector<GlycanIndices> glycositeIndices;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            size_t site = codeUtils::indexOf(residues, glycosites[n].GetResidue());
            size_t moleculeIndex =
                codeUtils::indexOf(molecules, codeUtils::erratic_cast<cds::Molecule*>(glycosites[n].GetGlycan()));
            auto& linkages                 = glycosidicLinkages[n];
            std::vector<size_t> linkageIds = codeUtils::indexVectorWithOffset(residueLinkages.size(), linkages);
            for (auto& linkage : linkages)
            {
                auto& linkageDihedrals = linkage.rotatableDihedrals;
                std::vector<size_t> dihedralIndices =
                    codeUtils::indexVectorWithOffset(rotatableDihedralIndices.size(), linkageDihedrals);
                codeUtils::insertInto(rotatableDihedralCurrentShape,
                                      cds::currentShape(linkage.rotatableDihedrals, linkage.dihedralMetadata));
                for (auto& dihedral : linkageDihedrals)
                {
                    std::array<size_t, 4> dihedralAtoms;
                    for (size_t n = 0; n < 4; n++)
                    {
                        dihedralAtoms[n] = codeUtils::indexOf(atoms, dihedral.atoms[n]);
                    }
                    rotatableDihedralIndices.push_back({dihedralAtoms, indexOfAtoms(dihedral.movingAtoms)});
                }
                auto onlyThisMolecule = [&](const std::vector<size_t>& indices)
                {
                    std::vector<size_t> result;
                    result.reserve(indices.size());
                    for (size_t index : indices)
                    {
                        if (graphIndices.residueMolecule[index] == moleculeIndex)
                        {
                            result.push_back(index);
                        }
                    }
                    return result;
                };
                std::vector<size_t> nonReducing = onlyThisMolecule(indexOfResidues(linkage.nonReducingOverlapResidues));
                std::vector<size_t> reducing    = onlyThisMolecule(indexOfResidues(linkage.reducingOverlapResidues));
                size_t firstResidue             = codeUtils::indexOf(residues, linkage.link.residues.first);
                size_t secondResidue            = codeUtils::indexOf(residues, linkage.link.residues.second);
                const std::vector<size_t>& adjacencies = residueGraph.nodes.nodeAdjacencies[firstResidue];
                size_t edgeN                           = codeUtils::indexOf(adjacencies, secondResidue);
                if (edgeN >= adjacencies.size())
                {
                    throw std::runtime_error("no residue adjacency");
                }
                size_t edgeId                    = residueGraph.nodes.edgeAdjacencies[firstResidue][edgeN];
                std::array<size_t, 2> residueIds = residueGraph.edges.nodeAdjacencies[edgeId];
                std::array<size_t, 2> atomIds    = atomGraph.edges.nodeAdjacencies[residueGraph.edges.indices[edgeId]];
                bool direction                   = graphIndices.atomResidue[atomIds[0]] == residueIds[1];
                std::array<std::vector<bool>, 2> bondedAtoms = {closelyBondedAtoms(residueIds[0], atomIds[direction]),
                                                                closelyBondedAtoms(residueIds[1], atomIds[!direction])};

                residueLinkages.push_back({edgeId, dihedralIndices, bondedAtoms, nonReducing, reducing});
                linkageRotamerTypes.push_back(linkage.rotamerType);
                linkageMetadata.push_back(linkage.dihedralMetadata);
                linkageOverlapBonds.push_back({
                    {residueIds, bondedAtoms}
                });
                linkageBranching.push_back(linkage.rotatableDihedrals[0].isBranchingLinkage);
            }
            moleculeTypes[moleculeIndex] = MoleculeType::glycan;
            glycositeIndices.push_back({site, moleculeIndex, linkageIds});
        }
        std::vector<double> residueOverlapWeight;
        residueOverlapWeight.reserve(residues.size());
        for (size_t molecule : graphIndices.residueMolecule)
        {
            residueOverlapWeight.push_back(moleculeTypes[molecule] == MoleculeType::protein ? overlapWeight.protein
                                                                                            : overlapWeight.glycan);
        }

        std::vector<cds::Sphere> residueBoundingSpheres =
            boundingSpheresOf(atomBoundingSpheres, residueGraph.nodes.elements);
        ResidueData residueData {cds::residueNames(residues), cds::residueTypes(residues),
                                 cds::residueNumbers(residues), residueOverlapWeight, residueBoundingSpheres};
        std::vector<cds::Sphere> moleculeBounds =
            boundingSpheresOf(residueBoundingSpheres, moleculeGraph.nodes.elements);

        MoleculeData moleculeData {moleculeTypes, moleculeBounds};
        RotatableDihedralData rotatableDihedralData {rotatableDihedralCurrentShape};
        ResidueLinkagedata residueLinkageData {linkageRotamerTypes, linkageMetadata, linkageOverlapBonds,
                                               linkageBranching};
        GlycanData glycanData {std::vector<bool>(glycosites.size(), true)};

        AssemblyData data {atomData, residueData, moleculeData, rotatableDihedralData, residueLinkageData, glycanData};

        std::vector<size_t> proteinMolecules;
        for (size_t n = 0; n < moleculeTypes.size(); n++)
        {
            if (moleculeTypes[n] == MoleculeType::protein)
            {
                proteinMolecules.push_back(n);
            }
        }

        AssemblyGraphs graphs {graphIndices,    atomGraph,        residueGraph,
                               moleculeGraph,   proteinMolecules, rotatableDihedralIndices,
                               residueLinkages, glycositeIndices};
        return GlycoproteinAssembly {graphs, data};
    }
} // namespace glycoproteinBuilder
