#include "include/readers/Pdb/amberMdPrep.hpp"

#include "include/readers/Pdb/bondByDistance.hpp"
#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbFunctions.hpp"
#include "include/readers/Pdb/pdbResidue.hpp"

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
                    ppInfo.nonNaturalProteinResidues_.emplace_back(pdbResidueId(data, residueId));
                    return true;
                }
            }
            return false;
        }
    } // namespace pdb
} // namespace gmml
