#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "include/CentralDataStructure/FileFormats/offFileData.hpp"
#include "include/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"

namespace gmml
{
    namespace gpbuilder
    {
        OffFileData toOffFileData(
            const assembly::Graph& graph, const AssemblyData& data, const std::vector<Coordinate>& atomCoordinates);

        PdbFileData toPdbFileData(
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
