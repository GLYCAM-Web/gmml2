#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/proteinBonding.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <vector>

namespace
{
    void setBondingForAminoAcid(cds::Residue* proteinRes)
    {
        for (auto& bondPair : biology::getBackboneBonding())
        {
            cds::Atom* firstAtom  = proteinRes->FindAtom(bondPair.first);
            cds::Atom* secondAtom = proteinRes->FindAtom(bondPair.second);
            if (firstAtom != nullptr && secondAtom != nullptr)
            {
                addBond(firstAtom, secondAtom);
            }
        }
        for (auto& bondPair : biology::getSidechainBonding(proteinRes->getName()))
        {
            cds::Atom* firstAtom  = proteinRes->FindAtom(bondPair.first);
            cds::Atom* secondAtom = proteinRes->FindAtom(bondPair.second);
            if (firstAtom != nullptr && secondAtom != nullptr)
            {
                addBond(firstAtom, secondAtom);
            }
        }
        return;
    }

    bool autoConnectSuccessiveResidues(cds::Residue* cTermRes, cds::Residue* nTermRes)
    {
        cds::Atom* cAtom = cTermRes->FindAtom("C");
        cds::Atom* nAtom = nTermRes->FindAtom("N");
        if (cAtom != nullptr && nAtom != nullptr && isWithinBondingDistance(cAtom, nAtom))
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
