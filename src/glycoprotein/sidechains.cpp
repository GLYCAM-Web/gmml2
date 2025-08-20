#include "include/glycoprotein/sidechains.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/matrix.hpp"
#include "include/geometry/orientation.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycanShape.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/graph/graphFunctions.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/metadata/sidechainRotamers.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
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

            std::array<Coordinate, 4> sidechainDihedralCoordinates(
                const std::vector<Sphere>& coords, const std::array<size_t, 4>& indices)
            {
                auto coord = [&](size_t n) { return coords[n].center; };
                return {coord(indices[0]), coord(indices[1]), coord(indices[2]), coord(indices[3])};
            }

            std::vector<size_t> sidechainMovingAtoms(const AssemblyData& data, size_t residue)
            {
                // perhaps this should be its own variable rather than 0-th element lookup
                return data.residues.sidechainDihedrals[residue][0].movingAtoms;
            }

            void moveSphereCenters(
                std::vector<Sphere>& coords, const Matrix4x4& matrix, const std::vector<size_t>& toMove)
            {
                for (size_t n : toMove)
                {
                    coords[n].center = matrix * coords[n].center;
                }
            }

            void setSidechainRotation(
                std::vector<Sphere>& coords,
                const std::vector<SidechainDihedral>& dihedrals,
                const SidechainRotation& rotation)
            {
                for (size_t k = 0; k < dihedrals.size(); k++)
                {
                    Matrix4x4 matrix = rotationTo(
                        sidechainDihedralCoordinates(coords, dihedrals[k].atoms),
                        constants::toRadians(rotation.chi[k]));
                    moveSphereCenters(coords, matrix, dihedrals[k].movingAtoms);
                }
            }
        } // namespace

        bool sidechainHasGlycanOverlap(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<size_t>& glycans,
            size_t sidechainResidue)
        {
            double overlapTolerance = overlapSettings.tolerance;
            const std::vector<size_t>& sidechainAtoms = sidechainMovingAtoms(data, sidechainResidue);
            Sphere bounds = boundingSphere(util::indicesToValues(mutableData.bounds.atoms, sidechainAtoms));
            for (size_t glycanId : glycans)
            {
                size_t glycanMolecule = data.glycans.moleculeId[glycanId];
                if (spheresOverlap(overlapTolerance, bounds, mutableData.bounds.molecules[glycanMolecule]))
                {
                    for (size_t residue : moleculeResidues(graph, glycanMolecule))
                    {
                        if (spheresOverlap(overlapTolerance, bounds, mutableData.bounds.residues[residue]))
                        {
                            for (size_t atom : residueAtoms(graph, residue))
                            {
                                if (spheresOverlap(overlapTolerance, bounds, mutableData.bounds.atoms[atom]))
                                {
                                    for (size_t otherAtom : sidechainAtoms)
                                    {
                                        PotentialFactor factor = potentialFactor(
                                            overlapSettings.potentialTable,
                                            data.atoms.elements[atom],
                                            data.atoms.elements[otherAtom]);
                                        double overlap = overlapAmount(
                                            factor,
                                            overlapTolerance,
                                            mutableData.bounds.atoms[atom],
                                            mutableData.bounds.atoms[otherAtom]);
                                        if (compareOverlaps(overlap, overlapSettings.rejectionThreshold) == 1)
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

        void updateSidechainRotation(
            const SidechainRotamerData& sidechains,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t residue,
            size_t rotation)
        {
            setSidechainRotation(
                mutableData.bounds.atoms, data.residues.sidechainDihedrals[residue], sidechains.rotations[rotation]);
            mutableData.residueSidechainMoved[residue] = true;
            updateResidueBounds(graph, mutableData.bounds, residue);
            updateResidueMoleculeBounds(graph, mutableData.bounds, residue);
        }

        void restoreSidechainRotation(
            const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData, size_t residue)
        {
            for (size_t atom : sidechainMovingAtoms(data, residue))
            {
                mutableData.bounds.atoms[atom] = data.atoms.initialState[atom];
            }
            mutableData.residueSidechainMoved[residue] = false;
            updateResidueBounds(graph, mutableData.bounds, residue);
            updateResidueMoleculeBounds(graph, mutableData.bounds, residue);
        }

        void setSidechainToLowestOverlapState(
            const SidechainRotamerData& sidechains,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            const std::vector<size_t>& preference,
            size_t residue)
        {
            std::vector<size_t> potentialOverlaps = atomsWithinSidechainPotentialBounds(
                overlapSettings, graph, data, mutableData, data.atoms.includeInEachOverlapCheck, residue);
            restoreSidechainRotation(graph, data, mutableData, residue);
            double initialOverlap = overlapVectorSum(sidechainOverlap(
                overlapSettings,
                graph.source.indices,
                data,
                mutableData.bounds.atoms,
                sidechainMovingAtoms(data, residue),
                potentialOverlaps));
            IndexedOverlap bestRotation = lowestOverlapSidechainRotation(
                sidechains,
                overlapSettings,
                graph.source.indices,
                data,
                mutableData,
                preference,
                residue,
                potentialOverlaps);
            if (compareOverlaps(initialOverlap, bestRotation.overlap) > 0)
            {
                updateSidechainRotation(sidechains, graph, data, mutableData, residue, bestRotation.index);
            }
        }

        std::vector<double> sidechainOverlap(
            const OverlapSettings& overlapSettings,
            const assembly::Indices& indices,
            const AssemblyData& data,
            const std::vector<Sphere>& bounds,
            const std::vector<size_t>& atomsA,
            const std::vector<size_t>& atomsB)
        {
            std::vector<double> result(indices.atomCount, 0.0);
            for (size_t n : atomsA)
            {
                Element elementA = data.atoms.elements[n];
                for (size_t k : atomsB)
                {
                    Element elementB = data.atoms.elements[k];
                    PotentialFactor factor = potentialFactor(overlapSettings.potentialTable, elementA, elementB);
                    double overlap = overlapAmount(factor, overlapSettings.tolerance, bounds[n], bounds[k]);
                    result[indices.atomResidue[n]] += overlap;
                    result[indices.atomResidue[k]] += overlap;
                }
            }
            return result;
        }

        IndexedOverlap lowestOverlapSidechainRotation(
            const SidechainRotamerData& sidechains,
            const OverlapSettings& overlapSettings,
            const assembly::Indices& indices,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<size_t>& preference,
            size_t sidechainResidue,
            const std::vector<size_t>& otherAtoms)
        {
            const std::vector<size_t> rotations = data.residues.sidechainRotations[sidechainResidue];
            const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[sidechainResidue];
            const std::vector<size_t>& movingAtoms = sidechainMovingAtoms(data, sidechainResidue);
            std::vector<Sphere> coords = mutableData.bounds.atoms;
            std::vector<double> overlaps;
            overlaps.reserve(rotations.size());
            for (size_t n = 0; n < rotations.size(); n++)
            {
                coords = mutableData.bounds.atoms;
                setSidechainRotation(coords, dihedrals, sidechains.rotations[rotations[n]]);
                double overlap = overlapAboveThresholdSum(
                    overlapSettings.rejectionThreshold,
                    sidechainOverlap(overlapSettings, indices, data, coords, movingAtoms, otherAtoms));
                overlaps.push_back(overlap);
            }
            size_t lowestIndex = preference[0];
            for (size_t n = 1; n < preference.size(); n++)
            {
                size_t current = preference[n];
                if (compareOverlaps(overlaps[lowestIndex], overlaps[current]) > 0)
                {
                    lowestIndex = current;
                }
            }
            return IndexedOverlap {rotations[lowestIndex], overlaps[lowestIndex]};
        }

        std::vector<size_t> atomsWithinSidechainPotentialBounds(
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& mutableData,
            const std::vector<bool>& includedAtoms,
            size_t sidechainResidue)
        {
            const Sphere sidechainBounds = data.residues.sidechainPotentialBounds[sidechainResidue];
            std::vector<size_t> result;
            result.reserve(256); // some extra, don't need to be precise
            auto overlapsWithSidechainBounds = [&](const Sphere& sphere)
            { return spheresOverlap(overlapSettings.tolerance, sidechainBounds, sphere); };
            auto insertOverlappingAtomsOfMolecule = [&](size_t molecule)
            {
                if (overlapsWithSidechainBounds(mutableData.bounds.molecules[molecule]))
                {
                    for (size_t residue : moleculeResidues(graph, molecule))
                    {
                        if ((residue != sidechainResidue) &&
                            overlapsWithSidechainBounds(mutableData.bounds.residues[residue]))
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
            for (size_t molecule = 0; molecule < moleculeCount(graph.source); molecule++)
            {
                if (mutableData.moleculeIncluded[molecule])
                {
                    insertOverlappingAtomsOfMolecule(molecule);
                }
            }
            return result;
        }

        std::vector<std::vector<SidechainDihedral>> sidechainDihedrals(
            const AminoAcidTable& aminoAcidTable, const assembly::Graph& graph, const AssemblyData& data)
        {
            std::vector<std::vector<SidechainDihedral>> result(
                residueCount(graph.source), std::vector<SidechainDihedral> {});
            for (size_t n = 0; n < residueCount(graph.source); n++)
            {
                bool isProtein = data.residues.types[n] == ResidueType::Protein;
                bool hasExpectedAtoms = data.residues.hasAllExpectedAtoms[n];
                if (isProtein && hasExpectedAtoms)
                {
                    const std::string& residue = data.residues.names[n];
                    size_t aminoAcid = aminoAcidIndex(aminoAcidTable, residue);
                    const std::vector<std::array<std::string, 4>>& dihedrals =
                        aminoAcidTable.sidechainDihedralAtoms[aminoAcid];
                    if (dihedrals.size() > 0)
                    {
                        const std::vector<size_t>& atomIndices = residueAtoms(graph, n);
                        std::vector<std::string> atomNames = util::indicesToValues(data.atoms.names, atomIndices);
                        auto atom = [&](const std::string& name)
                        { return atomIndices[util::indexOf(atomNames, name)]; };

                        std::vector<SidechainDihedral> sidechainDihedrals;
                        sidechainDihedrals.reserve(dihedrals.size());
                        for (auto& dihedral : dihedrals)
                        {
                            std::function<size_t(size_t)> index = [&](size_t k) { return atom(dihedral[k]); };
                            std::array<size_t, 4> indices = dihedralIndices(index);
                            std::vector<bool> reachable = graph::reachableNodes(
                                graph.atoms,
                                util::indicesToBools(
                                    atomCount(graph.source), std::vector<size_t> {indices[1], indices[2]}),
                                indices[2]);
                            std::vector<size_t> movingAtoms = util::boolsToIndices(reachable);
                            sidechainDihedrals.push_back({indices, movingAtoms});
                        }
                        result[n] = sidechainDihedrals;
                    }
                }
            }
            return result;
        }

        SidechainRotationsAndWeights sidechainRotationsAndWeights(
            const assembly::Indices& indices, const AssemblyData& data, const SidechainRotamerData& sidechains)
        {
            std::vector<std::vector<size_t>> rotations(indices.residueCount, std::vector<size_t> {});
            std::vector<std::vector<double>> weights(indices.residueCount, std::vector<double> {});
            std::function<double(const size_t&)> weight = [&](const size_t& n)
            { return sidechains.rotations[n].probability; };
            for (size_t n = 0; n < indices.residueCount; n++)
            {
                if (!data.residues.sidechainDihedrals[n].empty())
                {
                    const std::string& residue = data.residues.names[n];
                    const std::string originalResidue = originalResidueName(residue);
                    double phi = constants::toDegrees(data.residues.phiAngles[n]);
                    double psi = constants::toDegrees(data.residues.psiAngles[n]);
                    std::vector<size_t> rotationIds = sidechainRotationIndices(sidechains, originalResidue, phi, psi);
                    rotations[n] = rotationIds;
                    weights[n] = util::vectorMap(weight, rotationIds);
                }
            }
            return {rotations, weights};
        }

        std::vector<Sphere> sidechainPotentialBounds(
            const assembly::Indices& indices,
            const AssemblyData& data,
            const MutableData& mutableData,
            const SidechainRotamerData& sidechains)
        {
            std::vector<Sphere> result(
                indices.residueCount,
                Sphere {
                    0.0, {0.0, 0.0, 0.0}
            });
            for (size_t n = 0; n < indices.residueCount; n++)
            {
                const std::vector<SidechainDihedral>& dihedrals = data.residues.sidechainDihedrals[n];
                const std::vector<size_t>& rotations = data.residues.sidechainRotations[n];
                if (!(dihedrals.empty() || rotations.empty()))
                {
                    const std::vector<size_t>& firstDihedralMovingAtoms = sidechainMovingAtoms(data, n);
                    Sphere bounds =
                        boundingSphere(util::indicesToValues(mutableData.bounds.atoms, firstDihedralMovingAtoms));
                    for (size_t rotationId : rotations)
                    {
                        std::vector<Sphere> coords = mutableData.bounds.atoms;
                        setSidechainRotation(coords, dihedrals, sidechains.rotations[rotationId]);
                        for (size_t atom : firstDihedralMovingAtoms)
                        {
                            bounds = boundingSphereIncluding(bounds, coords[atom]);
                        }
                    }
                    result[n] = bounds;
                }
            }
            return result;
        }

        std::vector<bool> partOfMovableSidechain(const assembly::Indices& indices, const AssemblyData& data)
        {
            std::vector<bool> result(indices.atomCount, false);
            for (size_t n = 0; n < indices.residueCount; n++)
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

        GlycoproteinAssembly addSidechainRotamers(
            const AminoAcidTable& aminoAcidTable, const SidechainRotamerData& sidechains, GlycoproteinAssembly assembly)
        {
            assembly.data.residues.sidechainDihedrals =
                sidechainDihedrals(aminoAcidTable, assembly.graph, assembly.data);
            SidechainRotationsAndWeights rotationsAndWeights =
                sidechainRotationsAndWeights(assembly.graph.source.indices, assembly.data, sidechains);
            assembly.data.residues.sidechainRotations = rotationsAndWeights.rotations;
            assembly.data.residues.sidechainRotationWeights = rotationsAndWeights.weights;
            assembly.data.residues.sidechainPotentialBounds = sidechainPotentialBounds(
                assembly.graph.source.indices, assembly.data, assembly.mutableData, sidechains);
            assembly.data.atoms.partOfMovableSidechain =
                partOfMovableSidechain(assembly.graph.source.indices, assembly.data);
            return assembly;
        }
    } // namespace gpbuilder
} // namespace gmml
