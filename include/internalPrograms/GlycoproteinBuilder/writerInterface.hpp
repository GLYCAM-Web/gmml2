#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/off/offFileData.hpp"
#include "include/pdb/pdbFileData.hpp"

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
