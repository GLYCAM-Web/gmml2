#ifndef INCLUDE_FILETYPE_PDB_PDBCHAIN_HPP
#define INCLUDE_FILETYPE_PDB_PDBCHAIN_HPP

#include "include/fileType/pdb/pdbData.hpp"

#include <sstream>
#include <string>

namespace gmml
{
    namespace pdb
    {
        void readChain(PdbData& data, size_t moleculeId, std::stringstream& stream_block);
        std::stringstream extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line);
        void tagTerminalResidues(PdbData& data, size_t moleculeId);
        size_t getNTerminal(const PdbData& data, size_t moleculeId);
        size_t getCTerminal(const PdbData& data, size_t moleculeId);
    } // namespace pdb
} // namespace gmml

#endif
