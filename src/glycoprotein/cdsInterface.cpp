#include "include/glycoprotein/cdsInterface.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/orientation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/casting.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        PartialLinkageData linkageData(
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const GraphIndexData& graphData,
            const assembly::Graph& graph,
            const std::vector<std::vector<ResidueLinkage>>& glycosidicLinkages)
        {
            const std::vector<Atom*>& atoms = graphData.objects.atoms;
            const std::vector<Residue*>& residues = graphData.objects.residues;

            auto indexOfAtoms = [&atoms](const std::vector<Atom*>& toFind)
            {
                std::vector<size_t> result;
                result.reserve(toFind.size());
                for (auto& res : toFind)
                {
                    result.push_back(util::indexOf(atoms, res));
                }
                return result;
            };

            std::vector<AngleWithMetadata> rotatableDihedralShape;
            std::vector<std::array<size_t, 4>> rotatableBondDihedralAtoms;
            std::vector<std::vector<size_t>> rotatableBondMovingAtoms;
            std::vector<size_t> residueLinkageEdgeId;
            std::vector<std::vector<size_t>> residueLinkageBonds;
            std::vector<RotamerType> linkageRotamerTypes;
            std::vector<std::vector<size_t>> dihedralMetadata;
            std::vector<size_t> currentDihedralMetadata;
            std::vector<bool> isGlycositeLinkage;
            std::vector<std::vector<size_t>> glycanLinkages;

            for (size_t n = 0; n < glycosidicLinkages.size(); n++)
            {
                const std::vector<ResidueLinkage>& linkages = glycosidicLinkages[n];
                std::vector<size_t> linkageIds = util::indexVectorWithOffset(residueLinkageEdgeId.size(), linkages);
                for (size_t k = 0; k < linkages.size(); k++)
                {
                    isGlycositeLinkage.push_back(k == 0);
                    const ResidueLinkage& linkage = linkages[k];
                    const std::vector<RotatableBond>& linkageRotatableBonds = linkage.rotatableBonds;
                    std::vector<size_t> bondIndices =
                        util::indexVectorWithOffset(rotatableBondDihedralAtoms.size(), linkageRotatableBonds);
                    util::insertInto(
                        rotatableDihedralShape,
                        currentShape(dihedralAngleDataTable, linkage.rotatableBonds, linkage.dihedralMetadata));
                    util::insertInto(dihedralMetadata, linkage.dihedralMetadata);
                    for (size_t q = 0; q < linkageRotatableBonds.size(); q++)
                    {
                        const RotatableBond& bond = linkageRotatableBonds[q];
                        std::array<size_t, 4> dihedralAtoms;
                        for (size_t i = 0; i < 4; i++)
                        {
                            dihedralAtoms[i] = util::indexOf(atoms, bond.dihedralAtoms[i]);
                        }
                        rotatableBondDihedralAtoms.push_back(dihedralAtoms);
                        rotatableBondMovingAtoms.push_back(indexOfAtoms(bond.movingAtoms));
                        currentDihedralMetadata.push_back(bond.currentMetadataIndex);
                    }
                    size_t firstResidue = util::indexOf(residues, linkage.link.residues.first);
                    size_t secondResidue = util::indexOf(residues, linkage.link.residues.second);
                    const std::vector<size_t>& adjacencies = graph.residues.nodes.nodeAdjacencies[firstResidue];
                    size_t edgeN = util::indexOf(adjacencies, secondResidue);
                    if (edgeN >= adjacencies.size())
                    {
                        throw std::runtime_error("no residue adjacency");
                    }
                    size_t edgeId = graph.residues.nodes.edgeAdjacencies[firstResidue][edgeN];
                    residueLinkageEdgeId.push_back(edgeId);
                    residueLinkageBonds.push_back(bondIndices);
                    linkageRotamerTypes.push_back(linkage.rotamerType);
                }
                glycanLinkages.push_back(linkageIds);
            }

            RotatableBondData rotatableBondData {
                rotatableBondDihedralAtoms, rotatableBondMovingAtoms, dihedralMetadata, rotatableDihedralShape};

            ResidueLinkageData residueLinkages {
                linkageRotamerTypes, isGlycositeLinkage, residueLinkageEdgeId, residueLinkageBonds};

            return {rotatableBondData, currentDihedralMetadata, residueLinkages, glycanLinkages};
        }

        GlycoproteinAssembly toGlycoproteinAssemblyStructs(
            const AminoAcidTable& aminoAcidTable,
            const util::SparseVector<double>& elementRadii,
            const pdb::PdbData& pdbData,
            const GraphIndexData& graphData,
            const assembly::Graph& graph,
            const std::vector<size_t>& glycanMoleculeIds,
            const std::vector<size_t>& glycositeIds,
            const PartialLinkageData& linkageData,
            bool excludeHydrogen)
        {
            const std::vector<Atom*>& atoms = graphData.objects.atoms;
            const std::vector<Residue*>& residues = graphData.objects.residues;

            std::vector<MoleculeType> moleculeTypes(graphData.indices.moleculeCount, MoleculeType::protein);
            util::setIndicesTo(
                moleculeTypes,
                glycanMoleculeIds,
                std::vector<MoleculeType>(glycanMoleculeIds.size(), MoleculeType::glycan));

            std::vector<std::string> atomNames = gmml::atomNames(atoms);
            std::vector<Sphere> atomBoundingSpheres =
                assembly::toAtomBounds(elementRadii, atomElements(atoms), atomCoordinates(atoms));
            assembly::Bounds bounds = assembly::toAssemblyBounds(graph, atomBoundingSpheres);
            std::vector<bool> partOfMovableSidechain(atoms.size(), false);
            std::vector<Element> atomElements = gmml::atomElements(atoms);
            std::function<bool(const size_t&)> nonHydrogen = [&](const size_t& n)
            { return atomElements[n] != Element::H; };
            std::vector<bool> allAtoms(atoms.size(), true);
            std::vector<bool> nonHydrogenAtoms = util::vectorMap(nonHydrogen, util::indexVector(atomElements));
            std::vector<bool> includeInOverlapCheck = excludeHydrogen ? nonHydrogenAtoms : allAtoms;

            std::vector<std::string> residueNames = gmml::residueNames(residues);
            std::vector<ResidueType> residueTypes = gmml::residueTypes(residues);
            std::vector<bool> residuesHaveAllExpectedAtoms(residues.size(), true);
            std::vector<double> phiAngles(residues.size(), 0.0);
            std::vector<double> psiAngles(residues.size(), 0.0);

            auto dihedralCoordinates = [&](const std::array<size_t, 4>& arr)
            {
                auto coord = [&](size_t n) { return atomBoundingSpheres[arr[n]].center; };
                return std::array<Coordinate, 4> {coord(0), coord(1), coord(2), coord(3)};
            };

            for (size_t n : util::indicesOfElement(residueTypes, ResidueType::Protein))
            {
                size_t aminoAcidIndex = gmml::aminoAcidIndex(aminoAcidTable, residueNames[n]);
                std::vector<size_t> nonHydrogenAtoms = util::vectorFilter(nonHydrogen, residueAtoms(graph, n));
                std::vector<std::string> nonHydrogenNames = util::indicesToValues(atomNames, nonHydrogenAtoms);
                auto atomIndex = [&](const std::string& str)
                { return nonHydrogenAtoms[util::indexOf(nonHydrogenNames, str)]; };
                size_t atomN = atomIndex("N");
                size_t atomCA = atomIndex("CA");
                size_t atomC = atomIndex("C");
                bool hasExpectedAtoms = (util::sorted(nonHydrogenNames) == aminoAcidTable.atomNames[aminoAcidIndex]);
                bool isNonTerminal = pdbData.residues.terminality[n] == pdb::NonTerminal;
                for (size_t residueBondId : graph.residues.nodes.edgeAdjacencies[n])
                {
                    bool direction = (n == graph.residues.edges.nodeAdjacencies[residueBondId][0]);
                    size_t atomBondId = residueEdgeToAtomEdgeIndex(graph, residueBondId);
                    const std::array<size_t, 2>& atomBond = graph.atoms.edges.nodeAdjacencies[atomBondId];
                    size_t thisAtom = atomBond[0 == direction];
                    size_t otherAtom = atomBond[1 == direction];
                    std::string otherName = atomNames[otherAtom];
                    if ((thisAtom == atomC) && (otherName == "N"))
                    {
                        psiAngles[n] = angle(dihedralCoordinates({atomN, atomCA, atomC, otherAtom}));
                    }
                    else if ((thisAtom == atomN) && (otherName == "C"))
                    {
                        phiAngles[n] = angle(dihedralCoordinates({otherAtom, atomN, atomCA, atomC}));
                    }
                }
                residuesHaveAllExpectedAtoms[n] = hasExpectedAtoms && isNonTerminal;
            }

            std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(residues.size());
            std::vector<std::vector<size_t>> sidechainRotations(residues.size());
            std::vector<std::vector<double>> sidechainWeights(residues.size());
            std::vector<Sphere> sidechainPotentialBounds(
                residues.size(),
                Sphere {
                    0.0, Coordinate {0.0, 0.0, 0.0}
            });

            auto serializeNonProtein = [](const std::vector<uint>& numbers, const std::vector<bool>& isProtein)
            {
                std::vector<uint> result = numbers;
                int maxNumber = util::vectorMax(uint(0), util::boolsToValues(numbers, isProtein));
                for (size_t n : util::boolsToIndices(util::vectorNot(isProtein)))
                {
                    maxNumber++;
                    result[n] = maxNumber;
                }
                return result;
            };

            std::function<bool(const MoleculeType&)> isProtein = [](const MoleculeType& type)
            { return type == MoleculeType::protein; };
            std::vector<MoleculeType> residueMoleculeTypes =
                util::indicesToValues(moleculeTypes, residueMolecules(graph.source.indices));
            std::vector<MoleculeType> atomMoleculeTypes =
                util::indicesToValues(residueMoleculeTypes, atomResidues(graph.source.indices));

            std::function<std::string(const size_t&)> chainId = [&](const size_t& n)
            {
                size_t residueId =
                    (residueMoleculeTypes[n] == MoleculeType::protein)
                        ? n
                        : glycositeIds[util::indexOf(glycanMoleculeIds, graphData.indices.residueMolecule[n])];
                return pdbData.residues.chainIds[residueId];
            };

            std::vector<std::string> chainIds = util::vectorMap(chainId, util::indexVector(residues));
            std::vector<bool> proteinResidue = util::vectorMap(isProtein, residueMoleculeTypes);
            std::vector<bool> proteinAtom = util::vectorMap(isProtein, atomMoleculeTypes);
            std::vector<uint> residueNumbers = serializeNonProtein(gmml::residueNumbers(residues), proteinResidue);
            std::vector<uint> atomNumbers = serializeNonProtein(gmml::atomNumbers(atoms), proteinAtom);

            AtomData atomData {
                atomNames,
                atomTypes(atoms),
                atomNumbers,
                util::serializedNumberVector(atoms.size()),
                atomAtomicNumbers(atoms),
                atomElementStrings(atoms),
                atomElements,
                atomCharges(atoms),
                atomBoundingSpheres,
                allAtoms,
                includeInOverlapCheck,
                includeInOverlapCheck,
                partOfMovableSidechain};

            ResidueData residueData {
                residueNames,
                residueTypes,
                residuesHaveAllExpectedAtoms,
                phiAngles,
                psiAngles,
                residueStringIds(residues),
                residueNumbers,
                util::serializedNumberVector(residues.size()),
                chainIds,
                sidechainDihedrals,
                sidechainRotations,
                sidechainWeights,
                sidechainPotentialBounds};

            MoleculeData moleculeData {moleculeTypes, util::vectorMap(isProtein, moleculeTypes)};
            GlycanData glycanData {glycositeIds, glycanMoleculeIds, linkageData.glycanLinkageIds};

            AssemblyData data {
                atomData,
                residueData,
                moleculeData,
                ResidueEdgeData {assembly::atomsCloseToResidueEdges(graph)},
                linkageData.bonds,
                linkageData.linkages,
                glycanData};

            std::vector<bool> moleculeIncluded(moleculeCount(graph.source), true);
            MutableData mutableData {
                bounds,
                linkageData.initialDihedralMetadata,
                moleculeIncluded,
                std::vector<bool>(residues.size(), false)};

            return GlycoproteinAssembly {graph, data, mutableData};
        }
    } // namespace gpbuilder
} // namespace gmml
