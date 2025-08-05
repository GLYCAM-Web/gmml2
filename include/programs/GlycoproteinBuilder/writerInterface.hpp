#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/programs/GlycoproteinBuilder/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        off::OffFileData toOffFileData(
            const assembly::Graph& graph, const AssemblyData& data, const std::vector<Coordinate>& atomCoordinates);

        pdb::PdbFileData toPdbFileData(
            const assembly::Indices& indices,
            const AssemblyData& data,
            const std::vector<Coordinate>& atomCoordinates,
            const std::vector<uint>& atomNumbers,
            const std::vector<uint>& residueNumbers,
            const std::vector<std::string>& chainIds,
            const std::vector<std::string>& headerLines);
    } // namespace gpbuilder
} // namespace gmml

#endif
