#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <vector>

namespace
{
    bool bondCloseAtoms(pdb::PdbData& data, size_t atom1, size_t atom2)
    {
        if (isWithinBondingDistance(data, atom1, atom2))
        {
            addBond(data, atom1, atom2);
            return true;
        }
        return false;
    }

    void bondCloseAtomsAndResidues(pdb::PdbData& data, const std::vector<std::vector<size_t>>& residueAtoms,
                                   size_t residue1, size_t residue2)
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
            std::string edgeName = pdb::residueStringId(data, residue1) + "--" + pdb::residueStringId(data, residue2);
            data.objects.residues[residue1]->addNeighbor(edgeName, data.objects.residues[residue2]);
        }
    }

    void bondResidueAtoms(pdb::PdbData& data, const std::vector<size_t>& residueIds,
                          const std::vector<std::vector<size_t>>& residueAtoms, size_t n)
    {
        size_t residue1                      = residueIds[n];
        const std::vector<size_t>& res1Atoms = residueAtoms[residue1];
        for (size_t k = n + 1; k < residueIds.size(); k++)
        {
            size_t residue2                      = residueIds[k];
            const std::vector<size_t>& res2Atoms = residueAtoms[residue2];
            if (cds::withinDistance(constants::residueDistanceOverlapCutoff, data.atoms.coordinates[res1Atoms[0]],
                                    data.atoms.coordinates[res2Atoms[0]]))
            {
                bondCloseAtomsAndResidues(data, residueAtoms, residue1, residue2);
            }
        }
    }

    std::vector<std::vector<size_t>> residueAtomVectors(const pdb::PdbData& data)
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

bool pdb::isWithinBondingDistance(const PdbData& data, size_t atomA, size_t atomB)
{
    double maxLength =
        MolecularMetadata::maxBondLengthByAtomType(data.atoms.elements[atomA], data.atoms.elements[atomB]);
    return withinDistance(maxLength, data.atoms.coordinates[atomA], data.atoms.coordinates[atomB]);
}

void pdb::bondAtomsByDistance(PdbData& data, const std::vector<size_t>& atoms)
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

void pdb::bondAtomsAndResiduesByDistance(PdbData& data)
{
    std::vector<std::vector<size_t>> residueAtoms = residueAtomVectors(data);
    std::vector<size_t> residueIds                = codeUtils::indexVector(residueAtoms);
    for (size_t n = 0; n < residueIds.size(); n++)
    { // First bond by distance for atoms within each residue
        bondAtomsByDistance(data, residueAtoms[residueIds[n]]);
        // Then for each residue, find other residues within reasonable residue distance.
        bondResidueAtoms(data, residueIds, residueAtoms, n);
    }
}

void pdb::distanceBondIntra(PdbData& data, const std::vector<size_t>& residues)
{
    std::vector<std::vector<size_t>> residueAtoms = codeUtils::indicesToValues(residueAtomVectors(data), residues);
    for (auto& atoms : residueAtoms)
    {
        bondAtomsByDistance(data, atoms);
    }
}

void pdb::distanceBondInter(PdbData& data, const std::vector<size_t>& residues)
{
    std::vector<std::vector<size_t>> residueAtoms = residueAtomVectors(data);
    for (size_t n = 0; n < residues.size(); n++)
    {
        bondResidueAtoms(data, residues, residueAtoms, n);
    }
}
