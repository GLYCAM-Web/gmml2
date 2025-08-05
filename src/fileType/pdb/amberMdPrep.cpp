#include "include/fileType/pdb/amberMdPrep.hpp"

#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"

#include <vector>

namespace gmml
{
    namespace pdb
    {
        bool checkForNonNaturalProteinResidues(
            const PdbData& data,
            const std::vector<size_t>& unknownResidues,
            size_t cAtom,
            PreprocessorInformation& ppInfo)
        {
            size_t atomCount = data.indices.atomCount;
            for (size_t residueId : unknownResidues)
            {
                size_t nAtom = findResidueAtom(data, residueId, "N");
                if (nAtom < atomCount && isWithinBondingDistance(data, nAtom, cAtom))
                {
                    ppInfo.nonNaturalProteinResidues.emplace_back(pdbResidueId(data, residueId));
                    return true;
                }
            }
            return false;
        }
    } // namespace pdb
} // namespace gmml
