#ifndef INCLUDE_PDB_PDBMODEL_HPP
#define INCLUDE_PDB_PDBMODEL_HPP

#include "include/pdb/pdbData.hpp"

#include <sstream>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        void readAssembly(PdbData& data, size_t assemblyId, Assembly& assembly, std::stringstream& stream_block);
        void Write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream);
        std::string extractChainId(const std::string& line);
        std::stringstream extractSingleChainFromRecordSection(
            std::stringstream& stream_block, std::string line, const std::string& initialChainID);
        void extractCoordinatesFromModel(PdbData& data, std::stringstream& stream_block, std::string line);
    } // namespace pdb
} // namespace gmml

#endif
