#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::ResidueLinkageShapePreference> randomLinkageShapePreference(
        const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId,
        std::function<std::vector<size_t>(GlycamMetadata::DihedralAngleDataVector metadataVector)> randomMetadata,
        std::function<double(GlycamMetadata::DihedralAngleData metadata)> randomAngle,
        bool freezeGlycositeResidueConformation);
} // namespace glycoproteinBuilder
#endif
