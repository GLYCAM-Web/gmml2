#ifndef INCLUDE_GLYCOPROTEIN_SHAPERANDOMIZATION_HPP
#define INCLUDE_GLYCOPROTEIN_SHAPERANDOMIZATION_HPP

#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        typedef std::function<std::vector<size_t>(pcg32&, const DihedralAngleDataTable&, const std::vector<size_t>&)>
            MetadataOrder;

        struct AngleSettings
        {
            double preferenceDeviation;
            double searchDeviation;
            double searchIncrement;
            uint halfIntervalSearches;
            MetadataOrder randomMetadata;
        };

        GlycanShapePreference randomLinkageShapePreference(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleTable,
            const AngleSettings& settings,
            const AssemblyData& data,
            const assembly::Bounds& bounds,
            size_t glycanId,
            std::function<double(pcg32&, const AngleSettings&, const DihedralAngleData& metadata)> randomAngle,
            bool freezeGlycositeResidueConformation);
    } // namespace gpbuilder
} // namespace gmml

#endif
