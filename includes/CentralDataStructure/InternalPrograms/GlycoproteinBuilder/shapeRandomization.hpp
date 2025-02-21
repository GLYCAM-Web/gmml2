#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::ResidueLinkageShapePreference> randomLinkageShapePreference(
        pcg32& rng, const AssemblyData& data, const assembly::Bounds& bounds, size_t glycanId,
        std::function<std::vector<size_t>(pcg32&, GlycamMetadata::DihedralAngleDataVector metadataVector)>
            randomMetadata,
        std::function<double(pcg32&, GlycamMetadata::DihedralAngleData metadata)> randomAngle,
        bool freezeGlycositeResidueConformation);
} // namespace glycoproteinBuilder
#endif
