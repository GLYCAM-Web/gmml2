#ifndef INCLUDE_PREPROCESS_PARAMETERMANAGER_HPP
#define INCLUDE_PREPROCESS_PARAMETERMANAGER_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/fileType/lib/libraryFile.hpp"
#include "include/fileType/pdb/pdbData.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace preprocess
    {
        struct ParameterManager
        {
            lib::LibraryData lib;
        };

        ParameterManager loadParameters(const std::string& baseDir);

        void createAtomsForResidue(
            const ParameterManager& parameters, Residue* queryResidue, const std::string glycamNameForResidue);
        void setAtomChargesForResidues(
            const ParameterManager& parameters, pdb::PdbData& data, const std::vector<size_t>& residueIds);
        bool setAtomChargesForResidue(const ParameterManager& parameters, pdb::PdbData& data, size_t residueId);
    } // namespace preprocess
} // namespace gmml

#endif
