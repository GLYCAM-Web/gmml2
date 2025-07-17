#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    cds::GlycanShapePreference randomLinkageShapePreference(
        pcg32& rng, const AngleSettings& settings, const AssemblyData& data, const assembly::Bounds& bounds,
        size_t glycanId,
        std::function<double(pcg32&, const AngleSettings&, const GlycamMetadata::DihedralAngleData& metadata)>
            randomAngle,
        bool freezeGlycositeResidueConformation);
} // namespace glycoproteinBuilder
#endif
