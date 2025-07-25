#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_SHAPERANDOMIZATION_HPP

#include "include/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "include/External_Libraries/PCG/pcg_random.h"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        GlycanShapePreference randomLinkageShapePreference(
            pcg32& rng,
            const AngleSettings& settings,
            const AssemblyData& data,
            const assembly::Bounds& bounds,
            size_t glycanId,
            std::function<double(pcg32&, const AngleSettings&, const DihedralAngleData& metadata)> randomAngle,
            bool freezeGlycositeResidueConformation);
    }
} // namespace gmml

#endif
