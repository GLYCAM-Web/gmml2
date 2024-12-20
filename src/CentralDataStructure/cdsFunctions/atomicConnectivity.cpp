#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

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

    void setBondingForAminoAcid(cds::Residue* proteinRes)
    {
        const MolecularMetadata::AminoAcid& aminoAcid = MolecularMetadata::aminoAcid(proteinRes->getName());
        addBonds(proteinRes, aminoAcid.zwitterion.bonds);
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

    void setProteinIntraConnectivity(std::vector<cds::Residue*> proteinResidues)
    {
        for (auto& aa : proteinResidues)
        {
            setBondingForAminoAcid(aa);
        }
        return;
    }

    void setProteinInterConnectivity(std::vector<cds::Residue*> proteinResidues)
    {
        if (proteinResidues.empty())
        {
            return;
        }
        setBondingForAminoAcid(proteinResidues.front()); // does the first one, and handles when size is 1.
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

std::vector<std::pair<size_t, size_t>> cds::atomPairVectorIndices(const std::vector<Atom*>& atoms,
                                                                  const std::vector<std::pair<Atom*, Atom*>>& pairs)
{
    std::vector<std::pair<size_t, size_t>> result;
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

void cds::setIntraConnectivity(std::vector<cds::Residue*> residues)
{
    setProteinIntraConnectivity(cdsSelections::selectResiduesByType(residues, ResidueType::Protein));
    bool invertSelection = true;
    cds::distanceBondIntra(cdsSelections::selectResiduesByType(residues, ResidueType::Protein, invertSelection));
}

void cds::setInterConnectivity(std::vector<cds::Residue*> residues)
{
    setProteinInterConnectivity(cdsSelections::selectResiduesByType(residues, ResidueType::Protein));
    bool invertSelection = true;
    cds::distanceBondInter(cdsSelections::selectResiduesByType(residues, ResidueType::Protein, invertSelection));
}
