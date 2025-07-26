#include "include/pdb/amberMdPrep.hpp"

#include "include/pdb/bondByDistance.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbFunctions.hpp"
#include "include/pdb/pdbResidue.hpp"

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
