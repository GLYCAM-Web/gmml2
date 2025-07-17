#ifndef INCLUDES_CENTRALDATASTRUCTURE_PARAMETERS_PARAMETERMANAGER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_PARAMETERS_PARAMETERMANAGER_HPP

#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"

#include <string>
#include <vector>

namespace cdsParameters
{
    struct ParameterManager
    {
        lib::LibraryData lib;
    };

    ParameterManager loadParameters(const std::string& baseDir);
    // Functions
    bool setAtomChargesForResidue(const ParameterManager& parameters, cds::Residue* queryResidue);
    void setAtomChargesForResidues(const ParameterManager& parameters, std::vector<cds::Residue*> queryResidues);
    void createAtomsForResidue(const ParameterManager& parameters, cds::Residue* queryResidue,
                               const std::string glycamNameForResidue);
    // Attributes
} // namespace cdsParameters
#endif
