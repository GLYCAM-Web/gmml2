#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_WRITERINTERFACE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"

namespace glycoproteinBuilder
{
    cds::OffFileData toOffFileData(const AssemblyGraphs& graphs, const AssemblyData& data);
    cds::PdbFileData toPdbFileData(const AssemblyGraphs& graphs, const AssemblyData& data);
} // namespace glycoproteinBuilder
#endif
