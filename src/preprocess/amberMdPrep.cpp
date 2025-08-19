#include "include/preprocess/amberMdPrep.hpp"

#include "include/assembly/assemblyIndices.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"

#include <vector>

namespace gmml
{
    namespace preprocess
    {
        bool checkForNonNaturalProteinResidues(
            const pdb::PdbData& data,
            const std::vector<size_t>& unknownResidues,
            size_t cAtom,
            PreprocessorInformation& ppInfo)
        {
            for (size_t residueId : unknownResidues)
            {
                size_t nAtom = findResidueAtom(data, residueId, "N");
                if (nAtom < atomCount(data.assembly) && isWithinBondingDistance(data, nAtom, cAtom))
                {
                    ppInfo.nonNaturalProteinResidues.emplace_back(pdbResidueId(data, residueId));
                    return true;
                }
            }
            return false;
        }
    } // namespace preprocess
} // namespace gmml
