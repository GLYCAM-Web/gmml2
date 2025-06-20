#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"

#include <string>
#include <sstream>
#include <vector>
#include <variant>

// #include <functional>
namespace pdb
{
    void readChain(PdbData& data, size_t moleculeId, std::stringstream& stream_block);
    std::stringstream extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line);
    void tagTerminalResidues(PdbData& data, size_t moleculeId);
    size_t getNTerminal(const PdbData& data, size_t moleculeId);
    size_t getCTerminal(const PdbData& data, size_t moleculeId);
} // namespace pdb
#endif
