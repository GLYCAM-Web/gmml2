#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include <vector>
#include <sstream>
#include <ostream>

namespace pdb
{
    void readAssembly(cds::Assembly& assembly, std::stringstream& stream_block);
    void ChangeResidueName(cds::Assembly& assembly, const std::string& selector, const std::string& newName);

    // Preprocessing functions
    void preProcessCysResidues(cds::Assembly& assembly, PreprocessorInformation& ppInfo);
    void preProcessHisResidues(cds::Assembly& assembly, PreprocessorInformation& ppInfo,
                               const PreprocessorOptions& inputOptions);
    void preProcessChainTerminals(cds::Assembly& assembly, PreprocessorInformation& ppInfo,
                                  const PreprocessorOptions& inputOptions);
    void preProcessGapsUsingDistance(cds::Assembly& assembly, PreprocessorInformation& ppInfo,
                                     const PreprocessorOptions& inputOptions);
    void preProcessMissingUnrecognized(cds::Assembly& assembly, PreprocessorInformation& ppInfo,
                                       const cdsParameters::ParameterManager& parmManager);
    //        void bondAtomsByDistance();
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(const cds::Assembly& assembly, std::ostream& out);
    void Write(const cds::Assembly& assembly, std::ostream& stream);

    std::string extractChainId(const std::string& line);
    std::stringstream extractSingleChainFromRecordSection(std::stringstream& stream_block, std::string line,
                                                          const std::string& initialChainID);
    void extractCoordinatesFromModel(cds::Assembly& assembly, std::stringstream& stream_block, std::string line);
} // namespace pdb
#endif
