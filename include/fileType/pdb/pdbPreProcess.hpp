#ifndef INCLUDE_FILETYPE_PDB_PDBPREPROCESS_HPP
#define INCLUDE_FILETYPE_PDB_PDBPREPROCESS_HPP

#include "include/carbohydrate/parameterManager.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbPreprocessorInputs.hpp"

#include <ostream>
#include <sstream>

namespace gmml
{
    namespace pdb
    {
        PreprocessorInformation preProcess(
            PdbFile& file, const ParameterManager& parameterManager, PreprocessorOptions options);
        void changeResidueName(
            PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName);
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
            PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const ParameterManager& parmManager);
    } // namespace pdb
} // namespace gmml

#endif
