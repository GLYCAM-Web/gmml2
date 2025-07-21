#ifndef INCLUDES_READERS_PDB_AMBERMDPREP_HPP
#define INCLUDES_READERS_PDB_AMBERMDPREP_HPP

#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbPreprocessorInputs.hpp"

#include <vector>

// The pdb preprocessing code that lives in pdbFile and pdbModel needs to become free functions that live in here.
// We will want to process things that aren't pdb files. Also the name has changed from preprocessing to AmberMdPrep as
// of ~May 2023
namespace gmml
{
    bool checkForNonNaturalProteinResidues(
        const pdb::PdbData& data,
        const std::vector<size_t>& unknownResidues,
        size_t cAtom,
        pdb::PreprocessorInformation& ppInfo);
}

#endif
