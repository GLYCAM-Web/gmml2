#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include <vector>
#include <sstream>
#include <ostream>

namespace pdb
{
    void readAssembly(PdbData& data, size_t assemblyId, cds::Assembly& assembly, std::stringstream& stream_block);
    void Write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream);
    std::string extractChainId(const std::string& line);
    std::stringstream extractSingleChainFromRecordSection(std::stringstream& stream_block, std::string line,
                                                          const std::string& initialChainID);
    void extractCoordinatesFromModel(PdbData& data, std::stringstream& stream_block, std::string line);
} // namespace pdb
#endif
