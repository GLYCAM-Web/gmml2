#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/sidechains.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphFunctions.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>
#include <string>
#include <array>
#include <vector>

namespace
{
    std::array<size_t, 4> dihedralIndices(std::function<size_t(size_t)>& index)
    {
        std::array<size_t, 4> result {0, 0, 0, 0};
        for (size_t k = 0; k < 4; k++)
        {
            result[k] = index(k);
        }
        return result;
    }

    std::array<cds::Coordinate, 4> sidechainDihedralCoordinates(const std::vector<cds::Sphere>& coords,
                                                                const std::array<size_t, 4>& indices)
    {
        auto coord = [&](size_t n)
        {
            return coords[n].center;
        };
        return {coord(indices[0]), coord(indices[1]), coord(indices[2]), coord(indices[3])};
    }

    void moveSphereCenters(std::vector<cds::Sphere>& coords, const cds::RotationMatrix& matrix,
                           const std::vector<size_t>& toMove)
    {
        for (size_t n : toMove)
        {
            coords[n].center = matrix * coords[n].center;
        }
    }

    void setSidechainRotation(std::vector<cds::Sphere>& coords,
                              const std::vector<glycoproteinBuilder::SidechainDihedral>& dihedrals,
                              const MolecularMetadata::SidechainRotation& rotation)
    {
        for (size_t k = 0; k < dihedrals.size(); k++)
        {
            cds::RotationMatrix matrix = cds::rotationTo(sidechainDihedralCoordinates(coords, dihedrals[k].atoms),
                                                         constants::toRadians(rotation.chi[k]));
            moveSphereCenters(coords, matrix, dihedrals[k].movingAtoms);
        }
    }
} // namespace

bool glycoproteinBuilder::sidechainHasGlycanOverlap(const assembly::Graph& graph, const AssemblyData& data,
                                                    const MutableData& mutableData, const std::vector<size_t>& glycans,
                                                    size_t sidechainResidue)
{
    const std::vector<size_t>& sidechainAtoms = data.residues.sidechainDihedrals[sidechainResidue][0].movingAtoms;
    cds::Sphere bounds = cds::boundingSphere(codeUtils::indicesToValues(mutableData.atomBounds, sidechainAtoms));
    for (size_t glycanId : glycans)
    {
        size_t glycanMolecule = data.glycans.moleculeId[glycanId];
        if (cds::spheresOverlap(constants::overlapTolerance, bounds, mutableData.moleculeBounds[glycanMolecule]))
        {
            for (size_t residue : moleculeResidues(graph, glycanMolecule))
            {
                if (cds::spheresOverlap(constants::overlapTolerance, bounds, mutableData.residueBounds[residue]))
                {
                    for (size_t atom : residueAtoms(graph, residue))
                    {
                        if (cds::spheresOverlap(constants::overlapTolerance, bounds, mutableData.atomBounds[atom]))
                        {
                            for (size_t otherAtom : sidechainAtoms)
                            {
                                if (cds::spheresOverlap(constants::overlapTolerance, mutableData.atomBounds[atom],
                                                        mutableData.atomBounds[otherAtom]))
                                {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

void glycoproteinBuilder::updateSidechainRotation(const MolecularMetadata::SidechainRotamerData& sidechains,
                                                  const assembly::Graph& graph, const AssemblyData& data,
                                                  MutableData& mutableData, size_t residue, size_t rotation)
{
    setSidechainRotation(mutableData.atomBounds, data.residues.sidechainDihedrals[residue],
                         sidechains.rotations[rotation]);
    updateResidueBounds(graph, mutableData, residue);
    updateResidueMoleculeBounds(graph, mutableData, residue);
}

void glycoproteinBuilder::restoreSidechainRotation(const assembly::Graph& graph, const AssemblyData& data,
                                                   MutableData& mutableData, size_t residue)
{
    for (size_t atom : data.residues.sidechainDihedrals[residue][0].movingAtoms)
    {
        mutableData.atomBounds[atom] = data.atoms.initialState[atom];
    }
    mutableData.residueSidechainMoved[residue] = false;
    updateResidueBounds(graph, mutableData, residue);
    updateResidueMoleculeBounds(graph, mutableData, residue);
}

void glycoproteinBuilder::setSidechainToLowestOverlapState(const MolecularMetadata::SidechainRotamerData& sidechains,
                                                           const assembly::Graph& graph, const AssemblyData& data,
                                                           MutableData& mutableData, size_t residue)
{
    std::vector<size_t> potentialOverlaps = atomsWithinSidechainPotentialBounds(graph, data, mutableData, residue);
    cds::Overlap initialOverlap =
        sidechainOverlap(graph, mutableData.atomBounds, data.residues.overlapWeights,
                         data.residues.sidechainDihedrals[residue][0].movingAtoms, potentialOverlaps);
    IndexedOverlap bestRotation =
        lowestOverlapSidechainRotation(sidechains, graph, data, mutableData, residue, potentialOverlaps);
    if (cds::compareOverlaps(initialOverlap, bestRotation.overlap) > 0)
    {
        updateSidechainRotation(sidechains, graph, data, mutableData, residue, bestRotation.index);
        mutableData.residueSidechainMoved[residue] = true;
    }
}

cds::Overlap glycoproteinBuilder::sidechainOverlap(const assembly::Graph& graph, const std::vector<cds::Sphere>& bounds,
                                                   const std::vector<double>& residueOverlapWeight,
                                                   const std::vector<size_t>& atomsA, const std::vector<size_t>& atomsB)
{
    cds::OverlapProperties properties {constants::clashWeightBase, constants::overlapTolerance};
    cds::Overlap overlap {0.0, 0.0};
    for (size_t n : atomsA)
    {
        double weight = residueOverlapWeight[graph.atomResidue[n]];
        for (size_t k : atomsB)
        {
            overlap += cds::overlapAmount(properties, bounds[n], bounds[k]) * weight;
        }
    }
    return overlap;
}

bool glycoproteinBuilder::sidechainInitialStateHasOverlap(const assembly::Graph& graph, const AssemblyData& data,
                                                          const std::vector<cds::Coordinate>& initialCoords,
                                                          MutableData& mutableData, size_t sidechainResidue)
{
    const std::vector<size_t>& sidechainAtoms = data.residues.sidechainDihedrals[sidechainResidue][0].movingAtoms;
    std::vector<cds::Sphere> coords           = mutableData.atomBounds;
    for (size_t n : sidechainAtoms)
    {
        coords[n].center = initialCoords[n];
    }
    cds::Sphere bounds               = cds::boundingSphere(codeUtils::indicesToValues(coords, sidechainAtoms));
    auto overlapsWithSidechainBounds = [&](const cds::Sphere& sphere)
    {
        return cds::spheresOverlap(constants::overlapTolerance, bounds, sphere);
    };
    for (size_t molecule = 0; molecule < graph.moleculeCount; molecule++)
    {
        if (overlapsWithSidechainBounds(mutableData.moleculeBounds[molecule]))
        {
            for (size_t residue : moleculeResidues(graph, molecule))
            {
                if ((residue != sidechainResidue) && overlapsWithSidechainBounds(mutableData.residueBounds[residue]))
                {
                    for (size_t atom : residueAtoms(graph, residue))
                    {
                        if (overlapsWithSidechainBounds(mutableData.atomBounds[atom]))
                        {
                            for (size_t n : sidechainAtoms)
                            {
                                if (cds::spheresOverlap(constants::overlapTolerance, coords[atom], coords[n]))
                                {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

glycoproteinBuilder::IndexedOverlap glycoproteinBuilder::lowestOverlapSidechainRotation(
    const MolecularMetadata::SidechainRotamerData& sidechains, const assembly::Graph& graph, const AssemblyData& data,
    MutableData& mutableData, size_t sidechainResidue, const std::vector<size_t>& otherAtoms)
{
    const std::vector<size_t>& rotations            = data.residues.sidechainRotations[sidechainResidue];
    const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[sidechainResidue];
    const std::vector<size_t>& movingAtoms          = dihedrals[0].movingAtoms;
    std::vector<cds::Sphere> coords                 = mutableData.atomBounds;
    std::vector<cds::Overlap> overlaps;
    overlaps.reserve(rotations.size());
    for (size_t rotation : rotations)
    {
        coords = mutableData.atomBounds;
        setSidechainRotation(coords, dihedrals, sidechains.rotations[rotation]);
        cds::Overlap overlap = sidechainOverlap(graph, coords, data.residues.overlapWeights, movingAtoms, otherAtoms);
        overlaps.push_back(overlap);
    }
    size_t lowestIndex = 0;
    for (size_t n = 1; n < rotations.size(); n++)
    {
        if (cds::compareOverlaps(overlaps[lowestIndex], overlaps[n]) > 0)
        {
            lowestIndex = n;
        }
    }
    return IndexedOverlap {rotations[lowestIndex], overlaps[lowestIndex]};
}

std::vector<size_t> glycoproteinBuilder::atomsWithinSidechainPotentialBounds(const assembly::Graph& graph,
                                                                             const AssemblyData& data,
                                                                             const MutableData& mutableData,
                                                                             size_t sidechainResidue)
{
    const cds::Sphere sidechainBounds = data.residues.sidechainPotentialBounds[sidechainResidue];
    std::vector<size_t> result;
    result.reserve(256); // some extra, don't need to be precise
    auto overlapsWithSidechainBounds = [&](const cds::Sphere& sphere)
    {
        return cds::spheresOverlap(constants::overlapTolerance, sidechainBounds, sphere);
    };
    auto insertOverlappingAtomsOfMolecule = [&](size_t molecule)
    {
        if (overlapsWithSidechainBounds(mutableData.moleculeBounds[molecule]))
        {
            for (size_t residue : moleculeResidues(graph, molecule))
            {
                if ((residue != sidechainResidue) && overlapsWithSidechainBounds(mutableData.residueBounds[residue]))
                {
                    for (size_t atom : residueAtoms(graph, residue))
                    {
                        if (overlapsWithSidechainBounds(mutableData.atomBounds[atom]))
                        {
                            result.push_back(atom);
                        }
                    }
                }
            }
        }
    };
    for (size_t molecule = 0; molecule < graph.moleculeCount; molecule++)
    {
        if (mutableData.moleculeIncluded[molecule])
        {
            insertOverlappingAtomsOfMolecule(molecule);
        }
    }
    return result;
}

std::vector<std::vector<glycoproteinBuilder::SidechainDihedral>>
glycoproteinBuilder::sidechainDihedrals(const assembly::Graph& graph, const AssemblyData& data)
{
    std::vector<std::vector<SidechainDihedral>> result(graph.residueCount, std::vector<SidechainDihedral> {});
    const std::vector<std::string>& aminoAcids = MolecularMetadata::aminoAcidNames();
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        bool isProtein        = data.residues.types[n] == cds::ResidueType::Protein;
        bool hasExpectedAtoms = data.residues.hasAllExpectedAtoms[n];
        if (isProtein && hasExpectedAtoms)
        {
            const std::string& residue                               = data.residues.names[n];
            size_t aminoAcid                                         = codeUtils::indexOf(aminoAcids, residue);
            const std::vector<std::array<std::string, 4>>& dihedrals = MolecularMetadata::aminoAcidDihedrals(aminoAcid);
            if (dihedrals.size() > 0)
            {
                const std::vector<size_t>& atomIndices = residueAtoms(graph, n);
                std::vector<std::string> atomNames     = codeUtils::indicesToValues(data.atoms.names, atomIndices);
                auto atom                              = [&](const std::string& name)
                {
                    return atomIndices[codeUtils::indexOf(atomNames, name)];
                };

                std::vector<SidechainDihedral> sidechainDihedrals;
                sidechainDihedrals.reserve(dihedrals.size());
                for (auto& dihedral : dihedrals)
                {
                    std::function<size_t(size_t)> index = [&](size_t k)
                    {
                        return atom(dihedral[k]);
                    };
                    std::array<size_t, 4> indices   = dihedralIndices(index);
                    std::vector<size_t> movingAtoms = graph::reachableNodes(
                        graph.atoms,
                        codeUtils::indicesToBools(graph.atomCount, std::vector<size_t> {indices[1], indices[2]}),
                        indices[2]);
                    sidechainDihedrals.push_back({indices, movingAtoms});
                }
                result[n] = sidechainDihedrals;
            }
        }
    }
    return result;
}

std::vector<std::vector<size_t>>
glycoproteinBuilder::sidechainRotations(const assembly::Graph& graph, const AssemblyData& data,
                                        const MolecularMetadata::SidechainRotamerData& sidechains)
{
    std::vector<std::vector<size_t>> result(graph.residueCount, std::vector<size_t> {});
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        if (!data.residues.sidechainDihedrals[n].empty())
        {
            const std::string& residue        = data.residues.names[n];
            const std::string originalResidue = MolecularMetadata::originalResidueName(residue);
            double phi                        = constants::toDegrees(data.residues.phiAngles[n]);
            double psi                        = constants::toDegrees(data.residues.psiAngles[n]);
            std::vector<size_t> rotationIds =
                MolecularMetadata::sidechainRotationIndices(sidechains, originalResidue, phi, psi);
            result[n] = rotationIds;
        }
    }
    return result;
}

std::vector<cds::Sphere>
glycoproteinBuilder::sidechainPotentialBounds(const assembly::Graph& graph, const AssemblyData& data,
                                              const MutableData& mutableData,
                                              const MolecularMetadata::SidechainRotamerData& sidechains)
{
    std::vector<cds::Sphere> result(graph.residueCount, cds::Sphere {
                                                            0.0, {0.0, 0.0, 0.0}
    });
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[n];
        const std::vector<size_t>& rotations            = data.residues.sidechainRotations[n];
        if (!(dihedrals.empty() || rotations.empty()))
        {
            const std::vector<size_t>& firstDihedralMovingAtoms = dihedrals[0].movingAtoms;
            cds::Sphere bounds =
                cds::boundingSphere(codeUtils::indicesToValues(mutableData.atomBounds, firstDihedralMovingAtoms));
            for (size_t rotationId : rotations)
            {
                std::vector<cds::Sphere> coords = mutableData.atomBounds;
                setSidechainRotation(coords, dihedrals, sidechains.rotations[rotationId]);
                for (size_t atom : firstDihedralMovingAtoms)
                {
                    bounds = cds::boundingSphereIncluding(bounds, coords[atom]);
                }
            }
            result[n] = bounds;
        }
    }
    return result;
}

std::vector<bool> glycoproteinBuilder::partOfMovableSidechain(const assembly::Graph& graph, const AssemblyData& data)
{
    std::vector<bool> result(graph.atomCount, false);
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[n];
        if (dihedrals.size() > 0)
        {
            for (size_t k : dihedrals[0].movingAtoms)
            {
                result[k] = true;
            }
        }
    }
    return result;
}

glycoproteinBuilder::GlycoproteinAssembly
glycoproteinBuilder::addSidechainRotamers(const MolecularMetadata::SidechainRotamerData& sidechains,
                                          GlycoproteinAssembly assembly)
{
    assembly.data.residues.sidechainDihedrals = sidechainDihedrals(assembly.graph, assembly.data);
    assembly.data.residues.sidechainRotations = sidechainRotations(assembly.graph, assembly.data, sidechains);
    assembly.data.residues.sidechainPotentialBounds =
        sidechainPotentialBounds(assembly.graph, assembly.data, assembly.mutableData, sidechains);
    assembly.data.atoms.partOfMovableSidechain = partOfMovableSidechain(assembly.graph, assembly.data);
    assembly.data.atoms.alwaysIncluded         = codeUtils::vectorNot(assembly.data.atoms.partOfMovableSidechain);
    return assembly;
}
