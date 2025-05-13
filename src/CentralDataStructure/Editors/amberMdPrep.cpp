#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"

bool amberMdPrep::checkForNonNaturalProteinResidues(const pdb::PdbData& data,
                                                    std::vector<cds::Residue*> unknownResidues, const cds::Atom* cAtom,
                                                    pdb::PreprocessorInformation& ppInfo)
{
    for (auto& unknownResidue : unknownResidues)
    {
        size_t residueId = codeUtils::indexOf(data.objects.residues, unknownResidue);
        auto nAtom       = unknownResidue->FindAtom("N");
        if (nAtom != nullptr && isWithinBondingDistance(nAtom, cAtom))
        {
            ppInfo.nonNaturalProteinResidues_.emplace_back(pdb::pdbResidueId(data, residueId));
            return true;
        }
    }
    return false;
}
