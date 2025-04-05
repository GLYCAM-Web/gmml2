#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <array>
#include <vector>
#include <stdexcept>

namespace
{
    void addBonds(cds::Residue* proteinRes, const std::vector<std::pair<std::string, std::string>>& bonds)
    {
        for (auto& bondPair : bonds)
        {
            cds::Atom* firstAtom  = proteinRes->FindAtom(bondPair.first);
            cds::Atom* secondAtom = proteinRes->FindAtom(bondPair.second);
            if (firstAtom != nullptr && secondAtom != nullptr)
            {
                addBond(firstAtom, secondAtom);
            }
        }
    }

    void bondHydrogenAtomsToClosestNonHydrogen(std::vector<cds::Atom*>& hydrogenAtoms,
                                               std::vector<cds::Atom*>& nonHydrogenAtoms)
    {
        double hydrogenCutoff = MolecularMetadata::hydrogenCovalentBondMaxLength();
        for (cds::Atom* hydrogen : hydrogenAtoms)
        {
            cds::Coordinate coord = hydrogen->coordinate();
            size_t closest        = 0;
            double minDistance    = cds::distance(coord, nonHydrogenAtoms[0]->coordinate());
            for (size_t n = 1; n < nonHydrogenAtoms.size(); n++)
            {
                double distance = cds::distance(coord, nonHydrogenAtoms[n]->coordinate());
                if (distance < minDistance)
                {
                    if (minDistance < hydrogenCutoff)
                    {
                        gmml::log(__LINE__, __FILE__, gmml::WAR,
                                  "Hydrogen atom within " + std::to_string(minDistance) +
                                      " Å of multiple non-hydrogen atoms in residue");
                    }
                    minDistance = distance;
                    closest     = n;
                }
            }
            if (minDistance < hydrogenCutoff)
            {
                addBond(nonHydrogenAtoms[closest], hydrogen);
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::WAR,
                          "Hydrogen atom not within " + std::to_string(hydrogenCutoff) +
                              " Å of any non-hydrogen atom in residue");
            }
        }
    }

    void setBondingForAminoAcid(const MolecularMetadata::AminoAcidTable& table, cds::Residue* proteinRes)
    {
        size_t aminoAcid                         = MolecularMetadata::aminoAcidIndex(table, proteinRes->getName());
        std::vector<cds::Atom*> nonHydrogenAtoms = proteinRes->getNonHydrogenAtoms();
        addBonds(proteinRes, table.bonds[aminoAcid]);
        addBonds(proteinRes, MolecularMetadata::carboxylBonds());
        if (nonHydrogenAtoms.size() > 0)
        {
            std::vector<cds::Atom*> hydrogenAtoms = proteinRes->getHydrogenAtoms();
            bondHydrogenAtomsToClosestNonHydrogen(hydrogenAtoms, nonHydrogenAtoms);
        }
    }

    bool autoConnectSuccessiveResidues(cds::Residue* cTermRes, cds::Residue* nTermRes)
    {
        cds::Atom* cAtom   = cTermRes->FindAtom("C");
        cds::Atom* oxtAtom = cTermRes->FindAtom("OXT");
        cds::Atom* nAtom   = nTermRes->FindAtom("N");
        if ((cAtom != nullptr) && (nAtom != nullptr) && (oxtAtom == nullptr) && isWithinBondingDistance(cAtom, nAtom))
        {
            addBond(cAtom, nAtom);
            cTermRes->addNeighbor(nTermRes->getStringId() + "-" + cTermRes->getStringId(), nTermRes);
            return true;
        }
        return false;
    }

    void setProteinIntraConnectivity(const MolecularMetadata::AminoAcidTable& table,
                                     std::vector<cds::Residue*> proteinResidues)
    {
        for (auto& aa : proteinResidues)
        {
            setBondingForAminoAcid(table, aa);
        }
        return;
    }

    void setProteinInterConnectivity(const MolecularMetadata::AminoAcidTable& table,
                                     std::vector<cds::Residue*> proteinResidues)
    {
        if (proteinResidues.empty())
        {
            return;
        }
        setBondingForAminoAcid(table, proteinResidues.front()); // does the first one, and handles when size is 1.
        cds::Residue* previousRes = proteinResidues.front();
        for (std::vector<cds::Residue*>::iterator it = proteinResidues.begin() + 1; it != proteinResidues.end(); ++it)
        {
            if (!autoConnectSuccessiveResidues(previousRes, *it))
            { // Automatically bond the N and C atoms of successive residues
                gmml::log(__LINE__, __FILE__, gmml::WAR,
                          "Gap detected between " + previousRes->getStringId() + " and " + (*it)->getStringId());
            }
            previousRes = *it;
        }
        return;
    }
} // namespace

std::vector<std::pair<int, int>> cds::atomPairNumbers(const std::vector<std::pair<Atom*, Atom*>>& pairs)
{
    std::vector<std::pair<int, int>> result;
    result.reserve(pairs.size());
    for (auto& pair : pairs)
    {
        result.push_back({pair.first->getNumber(), pair.second->getNumber()});
    }
    return result;
}

std::vector<std::array<size_t, 2>> cds::atomPairVectorIndices(const std::vector<Atom*>& atoms,
                                                              const std::vector<std::pair<Atom*, Atom*>>& pairs)
{
    std::vector<std::array<size_t, 2>> result;
    result.reserve(pairs.size());
    for (auto& pair : pairs)
    {
        result.push_back({atomVectorIndex(atoms, pair.first), atomVectorIndex(atoms, pair.second)});
    }
    return result;
}

std::vector<std::pair<cds::Atom*, cds::Atom*>> cds::atomPairsConnectedToOtherResidues(std::vector<Atom*> atoms)
{
    std::vector<std::pair<Atom*, Atom*>> foundAtoms;
    for (auto& atom : atoms)
    { // only "child" neighbors or we find same pair twice
        for (auto& neighbor : atom->getChildren())
        { // check if neighbor is not one of the atoms in this residue.
            if (!codeUtils::contains(atoms, neighbor))
            {
                foundAtoms.push_back({atom, neighbor});
            }
        }
    }
    return foundAtoms;
}

std::vector<std::pair<cds::Atom*, cds::Atom*>> cds::atomPairsConnectedToOtherResidues(std::vector<Residue*> residues)
{
    std::vector<std::pair<Atom*, Atom*>> foundAtoms;
    for (auto& residue : residues)
    {
        auto newPairs = atomPairsConnectedToOtherResidues(residue->getAtoms());
        codeUtils::insertInto(foundAtoms, newPairs);
    }
    return foundAtoms;
}

std::vector<cds::Atom*> cds::atomsConnectedToOtherResidues(std::vector<Atom*> atoms)
{
    std::vector<Atom*> foundAtoms;
    for (auto& atom : atoms)
    {
        bool connected = false;
        for (auto& neighbor : atom->getNeighbors())
        { // check if neighbor is not one of the atoms in this residue.
            connected = connected || !codeUtils::contains(atoms, neighbor);
        }
        if (connected)
        {
            foundAtoms.push_back(atom);
        }
    }
    return foundAtoms;
}

void cds::setIntraConnectivity(const MolecularMetadata::AminoAcidTable& table, std::vector<cds::Residue*> residues)
{
    setProteinIntraConnectivity(table, cdsSelections::selectResiduesByType(residues, ResidueType::Protein));
    bool invertSelection = true;
    cds::distanceBondIntra(cdsSelections::selectResiduesByType(residues, ResidueType::Protein, invertSelection));
}

void cds::setInterConnectivity(const MolecularMetadata::AminoAcidTable& table, std::vector<cds::Residue*> residues)
{
    setProteinInterConnectivity(table, cdsSelections::selectResiduesByType(residues, ResidueType::Protein));
    bool invertSelection = true;
    cds::distanceBondInter(cdsSelections::selectResiduesByType(residues, ResidueType::Protein, invertSelection));
}
