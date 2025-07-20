#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBPREPROCESS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBPREPROCESS_HPP

#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"

#include <ostream>
#include <sstream>
#include <vector>

namespace pdb
{
    void changeResidueName(PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName);
    void preProcessCysResidues(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo);
    void preProcessHisResidues(
        PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions);
    void modifyTerminal(PdbData& data, size_t residueId, const std::string& type);
    void preProcessChainTerminals(
        PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions);
    void insertCap(PdbData& data, size_t moleculeId, size_t refResidueId, const std::string& type);
    void preProcessGapsUsingDistance(
        PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions);

    void preProcessMissingUnrecognized(
        PdbData& data,
        size_t assemblyId,
        PreprocessorInformation& ppInfo,
        const cdsParameters::ParameterManager& parmManager);
} // namespace pdb
#endif
