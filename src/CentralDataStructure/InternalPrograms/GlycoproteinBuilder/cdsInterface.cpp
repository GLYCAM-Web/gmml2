#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"

#include <array>
#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly toGlycoproteinAssemblyStructs(std::vector<cds::Molecule*>& molecules,
                                                       std::vector<GlycosylationSite>& glycosites,
                                                       const OverlapWeight overlapWeight)
    {
        using MolecularMetadata::Element;

        cds::GraphIndexData graphIndices = cds::toIndexData(molecules);
        std::vector<bool> includedAtoms(graphIndices.atoms.size(), true);
        assembly::Graph graph                = cds::createAssemblyGraph(graphIndices, includedAtoms);
        std::vector<cds::Atom*>& atoms       = graphIndices.atoms;
        std::vector<cds::Residue*>& residues = graphIndices.residues;

        auto boundingSpheresOf =
            [](const std::vector<cds::Sphere>& spheres, const std::vector<std::vector<size_t>>& indexVector)
        {
            std::vector<cds::Sphere> result;
            result.reserve(indexVector.size());
            for (auto& indices : indexVector)
            {
                result.push_back(cds::boundingSphere(codeUtils::indicesToValues(spheres, indices)));
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
            const std::vector<size_t>& elements = residueAtoms(graph, residueId);
            std::vector<bool> result(elements.size(), false);
            for (size_t n = 0; n < elements.size(); n++)
            {
                const std::vector<size_t>& adjacencies = graph.atoms.nodes.nodeAdjacencies[atomId];
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
        std::vector<cds::AngleWithMetadata> rotatableDihedralShape;
        std::vector<RotatableDihedralIndices> rotatableDihedralIndices;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycamMetadata::RotamerType> linkageRotamerTypes;
        std::vector<GlycamMetadata::DihedralAngleDataVector> dihedralMetadata;
        std::vector<std::vector<cds::BondedResidueOverlapInput>> linkageOverlapBonds;
        std::vector<bool> linkageBranching;
        std::vector<bool> isGlycositeLinkage;
        std::vector<GlycanIndices> glycositeIndices;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            size_t site = codeUtils::indexOf(residues, glycosites[n].GetResidue());
            size_t moleculeIndex =
                codeUtils::indexOf(molecules, codeUtils::erratic_cast<cds::Molecule*>(glycosites[n].GetGlycan()));
            const std::vector<cds::ResidueLinkage>& linkages = glycosidicLinkages[n];
            std::vector<size_t> linkageIds = codeUtils::indexVectorWithOffset(residueLinkages.size(), linkages);
            for (size_t k = 0; k < linkages.size(); k++)
            {
                isGlycositeLinkage.push_back(k == 0);
                const cds::ResidueLinkage& linkage                          = linkages[k];
                const std::vector<cds::RotatableDihedral>& linkageDihedrals = linkage.rotatableDihedrals;
                std::vector<size_t> dihedralIndices =
                    codeUtils::indexVectorWithOffset(rotatableDihedralIndices.size(), linkageDihedrals);
                codeUtils::insertInto(rotatableDihedralShape,
                                      cds::currentShape(linkage.rotatableDihedrals, linkage.dihedralMetadata));
                codeUtils::insertInto(dihedralMetadata, linkage.dihedralMetadata);
                for (size_t q = 0; q < linkageDihedrals.size(); q++)
                {
                    const cds::RotatableDihedral& dihedral = linkageDihedrals[q];
                    std::array<size_t, 4> dihedralAtoms;
                    for (size_t i = 0; i < 4; i++)
                    {
                        dihedralAtoms[i] = codeUtils::indexOf(atoms, dihedral.atoms[i]);
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
                const std::vector<size_t>& adjacencies = graph.residues.nodes.nodeAdjacencies[firstResidue];
                size_t edgeN                           = codeUtils::indexOf(adjacencies, secondResidue);
                if (edgeN >= adjacencies.size())
                {
                    throw std::runtime_error("no residue adjacency");
                }
                size_t edgeId                    = graph.residues.nodes.edgeAdjacencies[firstResidue][edgeN];
                std::array<size_t, 2> residueIds = graph.residues.edges.nodeAdjacencies[edgeId];
                std::array<size_t, 2> atomIds =
                    graph.atoms.edges.nodeAdjacencies[residueEdgeToAtomEdgeIndex(graph, edgeId)];
                bool direction                               = graphIndices.atomResidue[atomIds[0]] == residueIds[1];
                std::array<std::vector<bool>, 2> bondedAtoms = {closelyBondedAtoms(residueIds[0], atomIds[direction]),
                                                                closelyBondedAtoms(residueIds[1], atomIds[!direction])};

                residueLinkages.push_back({edgeId, dihedralIndices, bondedAtoms, nonReducing, reducing});
                linkageRotamerTypes.push_back(linkage.rotamerType);
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

        std::vector<std::string> atomNames           = cds::atomNames(atoms);
        std::vector<cds::Sphere> atomBoundingSpheres = cds::atomCoordinatesWithRadii(atoms);
        std::vector<bool> partOfMovableSidechain(atoms.size(), false);
        AtomData atomData {atomNames,
                           cds::atomTypes(atoms),
                           cds::atomNumbers(atoms),
                           cds::serializedNumberVector(atoms.size()),
                           cds::atomAtomicNumbers(atoms),
                           cds::atomElements(atoms),
                           cds::atomCharges(atoms),
                           atomBoundingSpheres,
                           partOfMovableSidechain,
                           std::vector<bool>(atoms.size(), true),
                           std::vector<bool>(atoms.size(), true)};

        std::vector<std::string> residueNames      = cds::residueNames(residues);
        std::vector<cds::ResidueType> residueTypes = cds::residueTypes(residues);
        std::vector<cds::Sphere> residueBoundingSpheres =
            boundingSpheresOf(atomBoundingSpheres, graph.residues.nodes.elements);
        std::vector<bool> residuesHaveAllExpectedAtoms(residues.size(), true);
        std::vector<double> phiAngles(residues.size(), 0.0);
        std::vector<double> psiAngles(residues.size(), 0.0);
        std::vector<Element> atomElementEnums          = cds::atomElementEnums(atoms);
        std::function<bool(const size_t&)> nonHydrogen = [&](const size_t& n)
        {
            return atomElementEnums[n] != Element::H;
        };

        auto dihedralCoordinates = [&](const std::array<size_t, 4>& arr)
        {
            auto coord = [&](size_t n)
            {
                return atomBoundingSpheres[arr[n]].center;
            };
            return std::array<cds::Coordinate, 4> {coord(0), coord(1), coord(2), coord(3)};
        };

        for (size_t n = 0; n < residues.size(); n++)
        {
            if (residueTypes[n] == cds::ResidueType::Protein)
            {
                const MolecularMetadata::AminoAcid& aminoAcid = MolecularMetadata::aminoAcid(residueNames[n]);
                std::vector<size_t> nonHydrogenAtoms          = codeUtils::filter(nonHydrogen, residueAtoms(graph, n));
                std::vector<std::string> nonHydrogenNames     = codeUtils::indicesToValues(atomNames, nonHydrogenAtoms);
                auto atomIndex                                = [&](const std::string& str)
                {
                    return nonHydrogenAtoms[codeUtils::indexOf(nonHydrogenNames, str)];
                };
                size_t atomN          = atomIndex("N");
                size_t atomCA         = atomIndex("CA");
                size_t atomC          = atomIndex("C");
                bool hasExpectedAtoms = (codeUtils::sorted(nonHydrogenNames) == aminoAcid.atomNames);
                bool isNTerminal      = true;
                bool isCTerminal      = true;
                for (size_t residueBondId : graph.residues.nodes.edgeAdjacencies[n])
                {
                    bool direction    = (n == graph.residues.edges.nodeAdjacencies[residueBondId][0]);
                    size_t atomBondId = residueEdgeToAtomEdgeIndex(graph, residueBondId);
                    const std::array<size_t, 2>& atomBond = graph.atoms.edges.nodeAdjacencies[atomBondId];
                    size_t thisAtom                       = atomBond[0 == direction];
                    size_t otherAtom                      = atomBond[1 == direction];
                    std::string otherName                 = atomNames[otherAtom];
                    if ((thisAtom == atomC) && (otherName == "N"))
                    {
                        isCTerminal  = false;
                        psiAngles[n] = cds::angle(dihedralCoordinates({atomN, atomCA, atomC, otherAtom}));
                    }
                    else if ((thisAtom == atomN) && (otherName == "C"))
                    {
                        isNTerminal  = false;
                        phiAngles[n] = cds::angle(dihedralCoordinates({otherAtom, atomN, atomCA, atomC}));
                    }
                }
                residuesHaveAllExpectedAtoms[n] = hasExpectedAtoms && !(isNTerminal || isCTerminal);
            }
        }

        std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(residues.size());
        std::vector<std::vector<size_t>> sidechainRotations(residues.size());
        std::vector<cds::Sphere> sidechainPotentialBounds(residues.size(), cds::Sphere {
                                                                               0.0, cds::Coordinate {0.0, 0.0, 0.0}
        });

        ResidueData residueData {residueNames,
                                 residueTypes,
                                 residuesHaveAllExpectedAtoms,
                                 phiAngles,
                                 psiAngles,
                                 cds::residueStringIds(residues),
                                 cds::residueNumbers(residues),
                                 cds::serializedNumberVector(residues.size()),
                                 residueOverlapWeight,
                                 sidechainDihedrals,
                                 sidechainRotations,
                                 sidechainPotentialBounds};
        std::vector<cds::Sphere> moleculeBounds =
            boundingSpheresOf(residueBoundingSpheres, graph.molecules.nodes.elements);

        MoleculeData moleculeData {moleculeTypes};
        RotatableDihedralData rotatableDihedralData {dihedralMetadata, rotatableDihedralShape};
        ResidueLinkageData residueLinkageData {linkageRotamerTypes, linkageOverlapBonds, linkageBranching,
                                               isGlycositeLinkage};

        std::vector<size_t> proteinMolecules;
        for (size_t n = 0; n < moleculeTypes.size(); n++)
        {
            if (moleculeTypes[n] == MoleculeType::protein)
            {
                proteinMolecules.push_back(n);
            }
        }

        AssemblyIndices indices {proteinMolecules, rotatableDihedralIndices, residueLinkages, glycositeIndices};

        AssemblyData data {atomData, residueData, moleculeData, rotatableDihedralData, residueLinkageData, indices};

        std::vector<bool> glycanIncluded(glycosites.size(), true);
        MutableData mutableData {atomBoundingSpheres, residueBoundingSpheres, moleculeBounds, glycanIncluded,
                                 std::vector<bool>(residues.size(), false)};

        return GlycoproteinAssembly {graph, data, mutableData};
    }
} // namespace glycoproteinBuilder
