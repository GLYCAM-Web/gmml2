#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/elements.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly
    toGlycoproteinAssemblyStructs(const MolecularMetadata::AminoAcidTable& aminoAcidTable,
                                  const GlycamMetadata::DihedralAngleDataTable& dihedralAngleDataTable,
                                  const codeUtils::SparseVector<double>& elementRadii, const pdb::PdbData& pdbData,
                                  std::vector<cds::Molecule*>& molecules, std::vector<GlycosylationSite>& glycosites,
                                  std::vector<cdsCondensedSequence::Carbohydrate*>& glycans,
                                  const OverlapMultiplier overlapWeight, double overlapTolerance,
                                  double overlapRejectionThreshold, bool excludeHydrogen)
    {
        using MolecularMetadata::Element;

        cds::GraphIndexData graphIndices     = cds::toIndexData(molecules);
        assembly::Graph graph                = cds::createCompleteAssemblyGraph(graphIndices);
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

        std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
        for (auto& a : glycans)
        {
            glycosidicLinkages.push_back(a->GetGlycosidicLinkages());
        }

        std::vector<MoleculeType> moleculeTypes(graphIndices.molecules.size(), MoleculeType::protein);
        std::vector<cds::AngleWithMetadata> rotatableDihedralShape;
        std::vector<RotatableDihedralIndices> rotatableDihedralIndices;
        std::vector<ResidueLinkageIndices> residueLinkages;
        std::vector<GlycamMetadata::RotamerType> linkageRotamerTypes;
        std::vector<std::vector<size_t>> dihedralMetadata;
        std::vector<size_t> dihedralCurrentMetadata;
        std::vector<bool> isGlycositeLinkage;
        std::vector<size_t> glycanAttachmentResidue;
        std::vector<size_t> glycanMoleculeId;
        std::vector<std::vector<size_t>> glycanLinkages;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            size_t site          = codeUtils::indexOf(residues, glycosites[n].residue);
            size_t moleculeIndex = codeUtils::indexOf(molecules, codeUtils::erratic_cast<cds::Molecule*>(glycans[n]));
            const std::vector<cds::ResidueLinkage>& linkages = glycosidicLinkages[n];
            std::vector<size_t> linkageIds = codeUtils::indexVectorWithOffset(residueLinkages.size(), linkages);
            for (size_t k = 0; k < linkages.size(); k++)
            {
                isGlycositeLinkage.push_back(k == 0);
                const cds::ResidueLinkage& linkage                          = linkages[k];
                const std::vector<cds::RotatableDihedral>& linkageDihedrals = linkage.rotatableDihedrals;
                std::vector<size_t> dihedralIndices =
                    codeUtils::indexVectorWithOffset(rotatableDihedralIndices.size(), linkageDihedrals);
                codeUtils::insertInto(
                    rotatableDihedralShape,
                    cds::currentShape(dihedralAngleDataTable, linkage.rotatableDihedrals, linkage.dihedralMetadata));
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
                    dihedralCurrentMetadata.push_back(dihedral.currentMetadataIndex);
                }
                size_t firstResidue                    = codeUtils::indexOf(residues, linkage.link.residues.first);
                size_t secondResidue                   = codeUtils::indexOf(residues, linkage.link.residues.second);
                const std::vector<size_t>& adjacencies = graph.residues.nodes.nodeAdjacencies[firstResidue];
                size_t edgeN                           = codeUtils::indexOf(adjacencies, secondResidue);
                if (edgeN >= adjacencies.size())
                {
                    throw std::runtime_error("no residue adjacency");
                }
                size_t edgeId = graph.residues.nodes.edgeAdjacencies[firstResidue][edgeN];
                residueLinkages.push_back({edgeId, dihedralIndices});
                linkageRotamerTypes.push_back(linkage.rotamerType);
            }
            moleculeTypes[moleculeIndex] = MoleculeType::glycan;
            glycanAttachmentResidue.push_back(site);
            glycanMoleculeId.push_back(moleculeIndex);
            glycanLinkages.push_back(linkageIds);
        }

        std::vector<std::string> atomNames           = cds::atomNames(atoms);
        std::vector<cds::Sphere> atomBoundingSpheres = cds::atomCoordinatesWithRadii(elementRadii, atoms);
        std::vector<bool> partOfMovableSidechain(atoms.size(), false);
        std::vector<Element> atomElements = cds::atomElements(atoms);
        MolecularMetadata::validateElementsInPotentialTable(MolecularMetadata::potentialTable(), atomElements);
        std::function<bool(const size_t&)> nonHydrogen = [&](const size_t& n)
        {
            return atomElements[n] != Element::H;
        };
        std::vector<bool> allAtoms(atoms.size(), true);
        std::vector<bool> nonHydrogenAtoms = codeUtils::vectorMap(nonHydrogen, codeUtils::indexVector(atomElements));
        std::vector<bool> includeInOverlapCheck = excludeHydrogen ? nonHydrogenAtoms : allAtoms;

        std::vector<std::string> residueNames      = cds::residueNames(residues);
        std::vector<cds::ResidueType> residueTypes = cds::residueTypes(residues);
        std::vector<cds::Sphere> residueBoundingSpheres =
            boundingSpheresOf(atomBoundingSpheres, graph.residues.nodes.constituents);
        std::vector<bool> residuesHaveAllExpectedAtoms(residues.size(), true);
        std::vector<double> phiAngles(residues.size(), 0.0);
        std::vector<double> psiAngles(residues.size(), 0.0);

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
                size_t aminoAcidIndex = MolecularMetadata::aminoAcidIndex(aminoAcidTable, residueNames[n]);
                std::vector<size_t> nonHydrogenAtoms = codeUtils::vectorFilter(nonHydrogen, residueAtoms(graph, n));
                std::vector<std::string> nonHydrogenNames = codeUtils::indicesToValues(atomNames, nonHydrogenAtoms);
                auto atomIndex                            = [&](const std::string& str)
                {
                    return nonHydrogenAtoms[codeUtils::indexOf(nonHydrogenNames, str)];
                };
                size_t atomN  = atomIndex("N");
                size_t atomCA = atomIndex("CA");
                size_t atomC  = atomIndex("C");
                bool hasExpectedAtoms =
                    (codeUtils::sorted(nonHydrogenNames) == aminoAcidTable.atomNames[aminoAcidIndex]);
                bool isNTerminal = true;
                bool isCTerminal = true;
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
        std::vector<std::vector<double>> sidechainWeights(residues.size());
        std::vector<cds::Sphere> sidechainPotentialBounds(residues.size(), cds::Sphere {
                                                                               0.0, cds::Coordinate {0.0, 0.0, 0.0}
        });

        auto serializeNonProtein = [](const std::vector<int>& numbers, const std::vector<bool>& isProtein)
        {
            std::vector<int> result = numbers;
            int maxNumber           = codeUtils::vectorMax(0, codeUtils::boolsToValues(numbers, isProtein));
            for (size_t n : codeUtils::boolsToIndices(codeUtils::vectorNot(isProtein)))
            {
                maxNumber++;
                result[n] = maxNumber;
            }
            return result;
        };

        std::function<bool(const MoleculeType&)> isProtein = [](const MoleculeType& type)
        {
            return type == MoleculeType::protein;
        };
        std::vector<MoleculeType> residueMoleculeTypes =
            codeUtils::indicesToValues(moleculeTypes, graph.residueMolecule);
        std::vector<MoleculeType> atomMoleculeTypes =
            codeUtils::indicesToValues(residueMoleculeTypes, graph.atomResidue);

        std::function<std::string(const size_t&)> chainId = [&](const size_t& n)
        {
            size_t residueId =
                (residueMoleculeTypes[n] == MoleculeType::protein)
                    ? n
                    : glycanAttachmentResidue[codeUtils::indexOf(glycanMoleculeId, graphIndices.residueMolecule[n])];
            size_t pdbResidueId = codeUtils::indexOf(pdbData.indices.residues, residues[residueId]);
            return pdbData.residues.chainIds[pdbResidueId];
        };

        std::vector<std::string> chainIds = codeUtils::vectorMap(chainId, codeUtils::indexVector(residues));
        std::vector<bool> proteinResidue  = codeUtils::vectorMap(isProtein, residueMoleculeTypes);
        std::vector<bool> proteinAtom     = codeUtils::vectorMap(isProtein, atomMoleculeTypes);
        std::vector<int> residueNumbers   = serializeNonProtein(cds::residueNumbers(residues), proteinResidue);
        std::vector<int> atomNumbers      = serializeNonProtein(cds::atomNumbers(atoms), proteinAtom);

        AtomData atomData {atomNames,
                           cds::atomTypes(atoms),
                           atomNumbers,
                           cds::serializedNumberVector(atoms.size()),
                           cds::atomAtomicNumbers(atoms),
                           cds::atomElementStrings(atoms),
                           atomElements,
                           cds::atomCharges(atoms),
                           atomBoundingSpheres,
                           allAtoms,
                           includeInOverlapCheck,
                           includeInOverlapCheck,
                           partOfMovableSidechain};

        ResidueData residueData {residueNames,
                                 residueTypes,
                                 residuesHaveAllExpectedAtoms,
                                 phiAngles,
                                 psiAngles,
                                 cds::residueStringIds(residues),
                                 residueNumbers,
                                 cds::serializedNumberVector(residues.size()),
                                 chainIds,
                                 sidechainDihedrals,
                                 sidechainRotations,
                                 sidechainWeights,
                                 sidechainPotentialBounds};
        std::vector<cds::Sphere> moleculeBounds =
            boundingSpheresOf(residueBoundingSpheres, graph.molecules.nodes.constituents);

        MoleculeData moleculeData {moleculeTypes};
        GlycanData GlycanData {glycanAttachmentResidue, glycanMoleculeId, glycanLinkages};
        RotatableDihedralData rotatableDihedralData {dihedralMetadata, rotatableDihedralShape};
        ResidueLinkageData residueLinkageData {linkageRotamerTypes, isGlycositeLinkage};

        std::vector<size_t> proteinMolecules;
        for (size_t n = 0; n < moleculeTypes.size(); n++)
        {
            if (moleculeTypes[n] == MoleculeType::protein)
            {
                proteinMolecules.push_back(n);
            }
        }

        AssemblyIndices indices {proteinMolecules, rotatableDihedralIndices, residueLinkages};

        cds::MoleculeOverlapWeight equalOverlapWeight {std::vector<double>(molecules.size(), 1.0),
                                                       std::vector<double>(molecules.size(), 1.0)};
        cds::MoleculeOverlapWeight defaultOverlapWeight;
        defaultOverlapWeight.within.reserve(molecules.size());
        defaultOverlapWeight.between.reserve(molecules.size());
        for (size_t molecule : graphIndices.residueMolecule)
        {
            bool isProtein = moleculeTypes[molecule] == MoleculeType::protein;
            defaultOverlapWeight.within.push_back(isProtein ? std::pow(overlapWeight.protein, 2.0)
                                                            : overlapWeight.self);
            defaultOverlapWeight.between.push_back(isProtein ? overlapWeight.protein : overlapWeight.glycan);
        }

        AssemblyData data {atomData,
                           residueData,
                           moleculeData,
                           ResidueEdgeData {assembly::atomsCloseToResidueEdges(graph)},
                           GlycanData,
                           rotatableDihedralData,
                           residueLinkageData,
                           indices,
                           dihedralAngleDataTable,
                           MolecularMetadata::potentialTable(),
                           defaultOverlapWeight,
                           equalOverlapWeight,
                           overlapTolerance,
                           overlapRejectionThreshold};

        std::vector<bool> moleculeIncluded(graph.moleculeCount, true);
        assembly::Bounds bounds {atomBoundingSpheres, residueBoundingSpheres, moleculeBounds};
        MutableData mutableData {bounds, dihedralCurrentMetadata, moleculeIncluded,
                                 std::vector<bool>(residues.size(), false)};

        return GlycoproteinAssembly {graph, data, mutableData};
    }
} // namespace glycoproteinBuilder
