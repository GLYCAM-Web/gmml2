#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"

bool amberMdPrep::checkForNonNaturalProteinResidues(const pdb::PdbData& data,
                                                    const std::vector<size_t>& unknownResidues, size_t cAtom,
                                                    pdb::PreprocessorInformation& ppInfo)
{
    size_t atomCount = data.indices.atomCount;
    for (size_t residueId : unknownResidues)
    {
        size_t nAtom = pdb::findResidueAtom(data, residueId, "N");
        if (nAtom < atomCount && pdb::isWithinBondingDistance(data, nAtom, cAtom))
        {
            ppInfo.nonNaturalProteinResidues_.emplace_back(pdb::pdbResidueId(data, residueId));
            return true;
        }
    }
    return false;
}
