#ifndef INCLUDE_PDB_AMBERMDPREP_HPP
#define INCLUDE_PDB_AMBERMDPREP_HPP

#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbPreprocessorInputs.hpp"

#include <vector>

// The pdb preprocessing code that lives in pdbFile and pdbModel needs to become free functions that live in here.
// We will want to process things that aren't pdb files. Also the name has changed from preprocessing to AmberMdPrep as
// of ~May 2023
namespace gmml
{
    namespace pdb
    {
        bool checkForNonNaturalProteinResidues(
            const PdbData& data,
            const std::vector<size_t>& unknownResidues,
            size_t cAtom,
            PreprocessorInformation& ppInfo);
    }
} // namespace gmml

#endif
