#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"

namespace glycoproteinBuilder
{
    cds::OffFileData toOffFileData(
        const assembly::Graph& graph, const AssemblyData& data, const std::vector<cds::Coordinate>& atomCoordinates);

    cds::PdbFileData toPdbFileData(
        const assembly::Indices& indices,
        const AssemblyData& data,
        const std::vector<cds::Coordinate>& atomCoordinates,
        const std::vector<uint>& atomNumbers,
        const std::vector<uint>& residueNumbers,
        const std::vector<std::string>& chainIds,
        const std::vector<std::string>& headerLines);
} // namespace glycoproteinBuilder
#endif
