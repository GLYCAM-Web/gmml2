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

    std::vector<size_t> sidechainMovingAtoms(const glycoproteinBuilder::AssemblyData& data, size_t residue)
    {
        // perhaps this should be its own variable rather than 0-th element lookup
        return data.residues.sidechainDihedrals[residue][0].movingAtoms;
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
    double overlapTolerance                   = data.overlapTolerance;
    const std::vector<size_t>& sidechainAtoms = sidechainMovingAtoms(data, sidechainResidue);
    cds::Sphere bounds = cds::boundingSphere(codeUtils::indicesToValues(mutableData.bounds.atoms, sidechainAtoms));
    for (size_t glycanId : glycans)
    {
        size_t glycanMolecule = data.glycans.moleculeId[glycanId];
        if (cds::spheresOverlap(overlapTolerance, bounds, mutableData.bounds.molecules[glycanMolecule]))
        {
            for (size_t residue : moleculeResidues(graph, glycanMolecule))
            {
                if (cds::spheresOverlap(overlapTolerance, bounds, mutableData.bounds.residues[residue]))
                {
                    for (size_t atom : residueAtoms(graph, residue))
                    {
                        if (cds::spheresOverlap(overlapTolerance, bounds, mutableData.bounds.atoms[atom]))
                        {
                            for (size_t otherAtom : sidechainAtoms)
                            {
                                if (cds::spheresOverlap(overlapTolerance, mutableData.bounds.atoms[atom],
                                                        mutableData.bounds.atoms[otherAtom]))
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
    setSidechainRotation(mutableData.bounds.atoms, data.residues.sidechainDihedrals[residue],
                         sidechains.rotations[rotation]);
    mutableData.residueSidechainMoved[residue] = true;
    updateResidueBounds(graph, mutableData.bounds, residue);
    updateResidueMoleculeBounds(graph, mutableData.bounds, residue);
}

void glycoproteinBuilder::restoreSidechainRotation(const assembly::Graph& graph, const AssemblyData& data,
                                                   MutableData& mutableData, size_t residue)
{
    for (size_t atom : sidechainMovingAtoms(data, residue))
    {
        mutableData.bounds.atoms[atom] = data.atoms.initialState[atom];
    }
    mutableData.residueSidechainMoved[residue] = false;
    updateResidueBounds(graph, mutableData.bounds, residue);
    updateResidueMoleculeBounds(graph, mutableData.bounds, residue);
}

void glycoproteinBuilder::setSidechainToLowestOverlapState(const MolecularMetadata::SidechainRotamerData& sidechains,
                                                           const assembly::Graph& graph, const AssemblyData& data,
                                                           MutableData& mutableData,
                                                           const std::vector<size_t>& preference, size_t residue)
{
    std::vector<size_t> potentialOverlaps =
        atomsWithinSidechainPotentialBounds(graph, data, mutableData, data.atoms.includeInEachOverlapCheck, residue);
    restoreSidechainRotation(graph, data, mutableData, residue);
    cds::Overlap initialOverlap = cds::overlapVectorSum(sidechainOverlap(
        graph, data, mutableData.bounds.atoms, sidechainMovingAtoms(data, residue), potentialOverlaps));
    IndexedOverlap bestRotation =
        lowestOverlapSidechainRotation(sidechains, graph, data, mutableData, preference, residue, potentialOverlaps);
    if (cds::compareOverlaps(initialOverlap, bestRotation.overlap) > 0)
    {
        updateSidechainRotation(sidechains, graph, data, mutableData, residue, bestRotation.index);
    }
}

std::vector<cds::Overlap> glycoproteinBuilder::sidechainOverlap(const assembly::Graph& graph, const AssemblyData& data,
                                                                const std::vector<cds::Sphere>& bounds,
                                                                const std::vector<size_t>& atomsA,
                                                                const std::vector<size_t>& atomsB)
{
    const cds::MoleculeOverlapWeight weight = data.defaultWeight;
    std::vector<cds::Overlap> result(graph.atomCount, {0.0, 0.0});
    for (size_t n : atomsA)
    {
        MolecularMetadata::Element elementA = data.atoms.elements[n];
        for (size_t k : atomsB)
        {
            double w = weight.between[graph.residueMolecule[graph.atomResidue[n]]] *
                       weight.between[graph.residueMolecule[graph.atomResidue[k]]];
            MolecularMetadata::Element elementB = data.atoms.elements[k];
            double scale                 = MolecularMetadata::potentialWeight(data.potentialTable, elementA, elementB);
            cds::Overlap overlap         = cds::overlapAmount(data.overlapTolerance, scale, bounds[n], bounds[k]) * w;
            result[graph.atomResidue[n]] += overlap;
            result[graph.atomResidue[k]] += overlap;
        }
    }
    return result;
}

glycoproteinBuilder::IndexedOverlap glycoproteinBuilder::lowestOverlapSidechainRotation(
    const MolecularMetadata::SidechainRotamerData& sidechains, const assembly::Graph& graph, const AssemblyData& data,
    const MutableData& mutableData, const std::vector<size_t>& preference, size_t sidechainResidue,
    const std::vector<size_t>& otherAtoms)
{
    const std::vector<size_t> rotations             = data.residues.sidechainRotations[sidechainResidue];
    const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[sidechainResidue];
    const std::vector<size_t>& movingAtoms          = sidechainMovingAtoms(data, sidechainResidue);
    std::vector<cds::Sphere> coords                 = mutableData.bounds.atoms;
    std::vector<cds::Overlap> overlaps;
    overlaps.reserve(rotations.size());
    for (size_t n = 0; n < rotations.size(); n++)
    {
        coords = mutableData.bounds.atoms;
        setSidechainRotation(coords, dihedrals, sidechains.rotations[rotations[n]]);
        cds::Overlap overlap = cds::overlapVectorSum(sidechainOverlap(graph, data, coords, movingAtoms, otherAtoms));
        overlaps.push_back(overlap);
    }
    size_t lowestIndex = preference[0];
    for (size_t n = 1; n < preference.size(); n++)
    {
        size_t current = preference[n];
        if (cds::compareOverlaps(overlaps[lowestIndex], overlaps[current]) > 0)
        {
            lowestIndex = current;
        }
    }
    return IndexedOverlap {rotations[lowestIndex], overlaps[lowestIndex]};
}

std::vector<size_t> glycoproteinBuilder::atomsWithinSidechainPotentialBounds(const assembly::Graph& graph,
                                                                             const AssemblyData& data,
                                                                             const MutableData& mutableData,
                                                                             const std::vector<bool>& includedAtoms,
                                                                             size_t sidechainResidue)
{
    const cds::Sphere sidechainBounds = data.residues.sidechainPotentialBounds[sidechainResidue];
    std::vector<size_t> result;
    result.reserve(256); // some extra, don't need to be precise
    auto overlapsWithSidechainBounds = [&](const cds::Sphere& sphere)
    {
        return cds::spheresOverlap(data.overlapTolerance, sidechainBounds, sphere);
    };
    auto insertOverlappingAtomsOfMolecule = [&](size_t molecule)
    {
        if (overlapsWithSidechainBounds(mutableData.bounds.molecules[molecule]))
        {
            for (size_t residue : moleculeResidues(graph, molecule))
            {
                if ((residue != sidechainResidue) && overlapsWithSidechainBounds(mutableData.bounds.residues[residue]))
                {
                    for (size_t atom : residueAtoms(graph, residue))
                    {
                        if (includedAtoms[atom] && overlapsWithSidechainBounds(mutableData.bounds.atoms[atom]))
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
glycoproteinBuilder::sidechainDihedrals(const MolecularMetadata::AminoAcidTable& aminoAcidTable,
                                        const assembly::Graph& graph, const AssemblyData& data)
{
    std::vector<std::vector<SidechainDihedral>> result(graph.residueCount, std::vector<SidechainDihedral> {});
    for (size_t n = 0; n < graph.residueCount; n++)
    {
        bool isProtein        = data.residues.types[n] == cds::ResidueType::Protein;
        bool hasExpectedAtoms = data.residues.hasAllExpectedAtoms[n];
        if (isProtein && hasExpectedAtoms)
        {
            const std::string& residue = data.residues.names[n];
            size_t aminoAcid           = MolecularMetadata::aminoAcidIndex(aminoAcidTable, residue);
            const std::vector<std::array<std::string, 4>>& dihedrals = aminoAcidTable.sidechainDihedralAtoms[aminoAcid];
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
                    std::array<size_t, 4> indices = dihedralIndices(index);
                    std::vector<bool> reachable   = graph::reachableNodes(
                        graph.atoms,
                        codeUtils::indicesToBools(graph.atomCount, std::vector<size_t> {indices[1], indices[2]}),
                        indices[2]);
                    std::vector<size_t> movingAtoms = codeUtils::boolsToIndices(reachable);
                    sidechainDihedrals.push_back({indices, movingAtoms});
                }
                result[n] = sidechainDihedrals;
            }
        }
    }
    return result;
}

glycoproteinBuilder::SidechainRotationsAndWeights
glycoproteinBuilder::sidechainRotationsAndWeights(const assembly::Graph& graph, const AssemblyData& data,
                                                  const MolecularMetadata::SidechainRotamerData& sidechains)
{
    std::vector<std::vector<size_t>> rotations(graph.residueCount, std::vector<size_t> {});
    std::vector<std::vector<double>> weights(graph.residueCount, std::vector<double> {});
    std::function<double(const size_t&)> weight = [&](const size_t& n)
    {
        return sidechains.rotations[n].probability;
    };
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
            rotations[n] = rotationIds;
            weights[n]   = codeUtils::vectorMap(weight, rotationIds);
        }
    }
    return {rotations, weights};
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
            const std::vector<size_t>& firstDihedralMovingAtoms = sidechainMovingAtoms(data, n);
            cds::Sphere bounds =
                cds::boundingSphere(codeUtils::indicesToValues(mutableData.bounds.atoms, firstDihedralMovingAtoms));
            for (size_t rotationId : rotations)
            {
                std::vector<cds::Sphere> coords = mutableData.bounds.atoms;
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
        if (data.residues.sidechainDihedrals[n].size() > 0)
        {
            for (size_t k : sidechainMovingAtoms(data, n))
            {
                result[k] = true;
            }
        }
    }
    return result;
}

glycoproteinBuilder::GlycoproteinAssembly
glycoproteinBuilder::addSidechainRotamers(const MolecularMetadata::AminoAcidTable& aminoAcidTable,
                                          const MolecularMetadata::SidechainRotamerData& sidechains,
                                          GlycoproteinAssembly assembly)
{
    assembly.data.residues.sidechainDihedrals = sidechainDihedrals(aminoAcidTable, assembly.graph, assembly.data);
    SidechainRotationsAndWeights rotationsAndWeights =
        sidechainRotationsAndWeights(assembly.graph, assembly.data, sidechains);
    assembly.data.residues.sidechainRotations       = rotationsAndWeights.rotations;
    assembly.data.residues.sidechainRotationWeights = rotationsAndWeights.weights;
    assembly.data.residues.sidechainPotentialBounds =
        sidechainPotentialBounds(assembly.graph, assembly.data, assembly.mutableData, sidechains);
    assembly.data.atoms.partOfMovableSidechain = partOfMovableSidechain(assembly.graph, assembly.data);
    return assembly;
}
