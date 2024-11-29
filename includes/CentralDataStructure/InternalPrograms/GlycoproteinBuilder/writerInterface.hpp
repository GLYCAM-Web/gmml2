#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"

namespace glycoproteinBuilder
{
    cds::OffFileData toOffFileData(const AssemblyGraphs& graphs, const AssemblyData& data,
                                   const std::vector<cds::Coordinate>& atomCoordinates);
    cds::PdbFileData toPdbFileData(const AssemblyGraphs& graphs, const AssemblyData& data,
                                   const std::vector<cds::Coordinate>& atomCoordinates,
                                   const std::vector<int>& atomNumbers, const std::vector<int>& residueNumbers);
} // namespace glycoproteinBuilder
#endif
