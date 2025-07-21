#include "include/readers/Pdb/bondByDistance.hpp"

#include "include/CentralDataStructure/residue.hpp"
#include "include/CentralDataStructure/residueTypes.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/metadata/elements.hpp"
#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbFunctions.hpp"
#include "include/readers/Pdb/pdbResidue.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        namespace
        {
            bool bondCloseAtoms(PdbData& data, size_t atom1, size_t atom2)
            {
                if (isWithinBondingDistance(data, atom1, atom2))
                {
                    addBond(data, atom1, atom2);
                    return true;
                }
                return false;
            }

            void bondCloseAtomsAndResidues(
                PdbData& data, const std::vector<std::vector<size_t>>& residueAtoms, size_t residue1, size_t residue2)
            {
                bool bonded = false;
                for (size_t atom1 : residueAtoms[residue1])
                {
                    for (size_t atom2 : residueAtoms[residue2])
                    {
                        bonded = bonded || bondCloseAtoms(data, atom1, atom2);
                    }
                }
                if (bonded)
                {
                    std::string edgeName = residueStringId(data, residue1) + "--" + residueStringId(data, residue2);
                    data.objects.residues[residue1]->addNeighbor(edgeName, data.objects.residues[residue2]);
                }
            }

            void bondResidueAtoms(
                PdbData& data,
                const std::vector<size_t>& residueIds,
                const std::vector<std::vector<size_t>>& residueAtoms,
                size_t n)
            {
                size_t residue1 = residueIds[n];
                const std::vector<size_t>& res1Atoms = residueAtoms[residue1];
                for (size_t k = n + 1; k < residueIds.size(); k++)
                {
                    size_t residue2 = residueIds[k];
                    const std::vector<size_t>& res2Atoms = residueAtoms[residue2];
                    if (withinDistance(
                            constants::residueDistanceOverlapCutoff,
                            data.atoms.coordinates[res1Atoms[0]],
                            data.atoms.coordinates[res2Atoms[0]]))
                    {
                        bondCloseAtomsAndResidues(data, residueAtoms, residue1, residue2);
                    }
                }
            }

            std::vector<std::vector<size_t>> residueAtomVectors(const PdbData& data)
            {
                std::vector<std::vector<size_t>> result(data.indices.residueCount, std::vector<size_t>());
                for (size_t n = 0; n < data.indices.atomCount; n++)
                {
                    if (data.indices.atomAlive[n])
                    {
                        result[data.indices.atomResidue[n]].push_back(n);
                    }
                }
                return result;
            }
        } // namespace

        bool isWithinBondingDistance(const PdbData& data, size_t atomA, size_t atomB)
        {
            double maxLength = maxBondLengthByAtomType(data.atoms.elements[atomA], data.atoms.elements[atomB]);
            return withinDistance(maxLength, data.atoms.coordinates[atomA], data.atoms.coordinates[atomB]);
        }

        void bondAtomsByDistance(PdbData& data, const std::vector<size_t>& atoms)
        {
            for (size_t n = 0; n < atoms.size(); n++)
            {
                for (size_t k = n + 1; k < atoms.size(); k++)
                {
                    size_t atom1 = atoms[n];
                    size_t atom2 = atoms[k];
                    if (isWithinBondingDistance(data, atom1, atom2))
                    {
                        addBond(data, atom1, atom2);
                    }
                }
            }
        }

        void bondAtomsAndResiduesByDistance(PdbData& data)
        {
            std::vector<std::vector<size_t>> residueAtoms = residueAtomVectors(data);
            std::vector<size_t> residueIds = util::indexVector(residueAtoms);
            for (size_t n = 0; n < residueIds.size(); n++)
            { // First bond by distance for atoms within each residue
                bondAtomsByDistance(data, residueAtoms[residueIds[n]]);
                // Then for each residue, find other residues within reasonable residue distance.
                bondResidueAtoms(data, residueIds, residueAtoms, n);
            }
        }

        void distanceBondIntra(PdbData& data, const std::vector<size_t>& residues)
        {
            std::vector<std::vector<size_t>> residueAtoms = util::indicesToValues(residueAtomVectors(data), residues);
            for (auto& atoms : residueAtoms)
            {
                bondAtomsByDistance(data, atoms);
            }
        }

        void distanceBondInter(PdbData& data, const std::vector<size_t>& residues)
        {
            std::vector<std::vector<size_t>> residueAtoms = residueAtomVectors(data);
            for (size_t n = 0; n < residues.size(); n++)
            {
                bondResidueAtoms(data, residues, residueAtoms, n);
            }
        }
    } // namespace pdb
} // namespace gmml
