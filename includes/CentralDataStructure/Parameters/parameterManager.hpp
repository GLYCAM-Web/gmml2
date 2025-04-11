#ifndef INCLUDES_CENTRALDATASTRUCTURE_PARAMETERS_PARAMETERMANAGER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_PARAMETERS_PARAMETERMANAGER_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <vector>
#include <unordered_map>

namespace cdsParameters
{
    struct ParameterManager
    {
        std::vector<std::string> residueNames;
        std::vector<cds::Residue*> residues;
        std::vector<cds::Molecule> libFiles;
    };

    ParameterManager loadParameters(const std::string& baseDir);
    // Functions
    cds::Residue* findParameterResidue(const ParameterManager& parameters, const std::string name);
    bool setAtomChargesForResidue(const ParameterManager& parameters, cds::Residue* queryResidue);
    void setAtomChargesForResidues(const ParameterManager& parameters, std::vector<cds::Residue*> queryResidues);
    void createAtomsForResidue(const ParameterManager& parameters, cdsCondensedSequence::ParsedResidue* queryResidue,
                               const std::string glycamNameForResidue);
    // Attributes
} // namespace cdsParameters
#endif
