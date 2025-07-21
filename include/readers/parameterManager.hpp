#ifndef INCLUDES_READERS_PARAMETERMANAGER_HPP
#define INCLUDES_READERS_PARAMETERMANAGER_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/readers/Lib/LibraryFile.hpp"

#include <string>
#include <vector>

namespace gmml
{
    struct ParameterManager
    {
        lib::LibraryData lib;
    };

    ParameterManager loadParameters(const std::string& baseDir);
    // Functions
    bool setAtomChargesForResidue(const ParameterManager& parameters, Residue* queryResidue);
    void setAtomChargesForResidues(const ParameterManager& parameters, std::vector<Residue*> queryResidues);
    void createAtomsForResidue(
        const ParameterManager& parameters, Residue* queryResidue, const std::string glycamNameForResidue);
    // Attributes
} // namespace gmml

#endif
