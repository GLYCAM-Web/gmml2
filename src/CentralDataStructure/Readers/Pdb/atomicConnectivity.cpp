#include "includes/CentralDataStructure/Readers/Pdb/atomicConnectivity.hpp"

#include "includes/Assembly/assemblyIndices.hpp"
#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"

#include <array>
#include <stdexcept>
#include <vector>

namespace
{
    void addBonds(pdb::PdbData& data, size_t proteinRes, const std::vector<std::pair<std::string, std::string>>& bonds)
    {
        size_t atomCount = data.indices.atomCount;
        for (auto& bondPair : bonds)
        {
            size_t firstAtom = pdb::findResidueAtom(data, proteinRes, bondPair.first);
            size_t secondAtom = pdb::findResidueAtom(data, proteinRes, bondPair.second);
            if (firstAtom < atomCount && secondAtom < atomCount)
            {
                pdb::addBond(data, firstAtom, secondAtom);
            }
        }
    }

    void bondHydrogenAtomsToClosestNonHydrogen(pdb::PdbData& data, size_t residueId)
    {
        std::vector<bool> isResidueAtom = assembly::isResidueAtom(data.indices, residueId);
        std::vector<bool> isHydrogen = codeUtils::vectorEquals(data.atoms.elements, MolecularMetadata::Element::H);
        std::vector<size_t> hydrogenAtoms =
            indicesOfLivingAtoms(data.indices, codeUtils::vectorAnd(isResidueAtom, isHydrogen));
        std::vector<size_t> nonHydrogenAtoms =
            indicesOfLivingAtoms(data.indices, codeUtils::vectorAnd(isResidueAtom, codeUtils::vectorNot(isHydrogen)));
        if (nonHydrogenAtoms.empty())
        {
            return;
        }
        double hydrogenCutoff = MolecularMetadata::hydrogenCovalentBondMaxLength();
        for (size_t hydrogen : hydrogenAtoms)
        {
            cds::Coordinate coord = data.atoms.coordinates[hydrogen];
            size_t closest = 0;
            double minDistance = cds::distance(coord, data.atoms.coordinates[nonHydrogenAtoms[0]]);
            for (size_t n = 1; n < nonHydrogenAtoms.size(); n++)
            {
                double distance = cds::distance(coord, data.atoms.coordinates[nonHydrogenAtoms[n]]);
                if (distance < minDistance)
                {
                    if (minDistance < hydrogenCutoff)
                    {
                        gmml::log(
                            __LINE__,
                            __FILE__,
                            gmml::WAR,
                            "Hydrogen atom within " + std::to_string(minDistance) +
                                " Å of multiple non-hydrogen atoms in residue");
                    }
                    minDistance = distance;
                    closest = n;
                }
            }
            if (minDistance < hydrogenCutoff)
            {
                addBond(data, nonHydrogenAtoms[closest], hydrogen);
            }
            else
            {
                gmml::log(
                    __LINE__,
                    __FILE__,
                    gmml::WAR,
                    "Hydrogen atom not within " + std::to_string(hydrogenCutoff) +
                        " Å of any non-hydrogen atom in residue");
            }
        }
    }

    void setBondingForAminoAcid(const MolecularMetadata::AminoAcidTable& table, pdb::PdbData& data, size_t residueId)
    {
        size_t aminoAcid = MolecularMetadata::aminoAcidIndex(table, data.residues.names[residueId]);
        addBonds(data, residueId, table.bonds[aminoAcid]);
        addBonds(data, residueId, MolecularMetadata::carboxylBonds());
        bondHydrogenAtomsToClosestNonHydrogen(data, residueId);
    }

    bool autoConnectSuccessiveResidues(pdb::PdbData& data, size_t cTermRes, size_t nTermRes)
    {
        size_t atomCount = data.indices.atomCount;
        size_t cAtom = pdb::findResidueAtom(data, cTermRes, "C");
        size_t oxtAtom = pdb::findResidueAtom(data, cTermRes, "OXT");
        size_t nAtom = pdb::findResidueAtom(data, nTermRes, "N");
        if ((cAtom < atomCount) && (nAtom < atomCount) && (oxtAtom == atomCount) &&
            pdb::isWithinBondingDistance(data, cAtom, nAtom))
        {
            addBond(data, cAtom, nAtom);
            return true;
        }
        return false;
    }

    void setProteinInterConnectivity(
        const MolecularMetadata::AminoAcidTable& table, pdb::PdbData& data, const std::vector<size_t>& residueIds)
    {
        if (residueIds.empty())
        {
            return;
        }
        setBondingForAminoAcid(table, data, residueIds[0]); // does the first one, and handles when size is 1.
        for (size_t n = 0; n < residueIds.size() - 1; n++)
        {
            size_t previousRes = residueIds[n];
            if (!autoConnectSuccessiveResidues(data, previousRes, residueIds[n + 1]))
            { // Automatically bond the N and C atoms of successive residues
                gmml::log(
                    __LINE__,
                    __FILE__,
                    gmml::WAR,
                    "Gap detected between " + pdb::residueStringId(data, previousRes) + " and " +
                        pdb::residueStringId(data, residueIds[n + 1]));
            }
        }
    }
} // namespace

void pdb::setIntraConnectivity(const MolecularMetadata::AminoAcidTable& table, PdbData& data)
{
    std::vector<bool> isProtein = codeUtils::vectorEquals(data.residues.types, cds::ResidueType::Protein);
    for (size_t residueId : codeUtils::boolsToIndices(isProtein))
    {
        setBondingForAminoAcid(table, data, residueId);
    }
    distanceBondIntra(data, codeUtils::boolsToIndices(codeUtils::vectorNot(isProtein)));
}

void pdb::setInterConnectivity(const MolecularMetadata::AminoAcidTable& table, PdbData& data)
{
    std::vector<bool> isProtein = codeUtils::vectorEquals(data.residues.types, cds::ResidueType::Protein);
    setProteinInterConnectivity(table, data, codeUtils::boolsToIndices(isProtein));
    distanceBondInter(data, codeUtils::boolsToIndices(codeUtils::vectorNot(isProtein)));
}
