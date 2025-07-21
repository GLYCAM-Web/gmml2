#include "include/readers/Pdb/amberMdPrep.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/readers/Pdb/bondByDistance.hpp"
#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbFunctions.hpp"
#include "include/readers/Pdb/pdbResidue.hpp"

namespace gmml
{
    bool checkForNonNaturalProteinResidues(
        const pdb::PdbData& data,
        const std::vector<size_t>& unknownResidues,
        size_t cAtom,
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
} // namespace gmml
