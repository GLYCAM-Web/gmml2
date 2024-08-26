#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"

bool amberMdPrep::checkForNonNaturalProteinResidues(std::vector<cds::Residue*> unknownResidues, const cds::Atom* cAtom,
                                                    pdb::PreprocessorInformation& ppInfo)
{
    for (auto& unknownResidue : unknownResidues)
    {
        auto nAtom = unknownResidue->FindAtom("N");
        if (nAtom != nullptr && isWithinBondingDistance(nAtom, cAtom))
        {
            ppInfo.nonNaturalProteinResidues_.emplace_back(unknownResidue->getId());
            return true;
        }
    }
    return false;
}
