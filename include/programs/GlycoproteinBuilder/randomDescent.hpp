#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/overlapCount.hpp"
#include "include/glycoprotein/shapeRandomization.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct GlycoproteinState
        {
            double overlap;
            std::vector<size_t> sitesWithOverlap;
            std::vector<size_t> sitesAboveOverlapThreshold;
            std::vector<GlycanShapePreference> preferences;
        };

        typedef std::function<void(
            const assembly::Graph&,
            const AssemblyData&,
            const assembly::Selection&,
            const AngleSettings&,
            const GlycanShapePreference&,
            MutableData&,
            size_t)>
            WiggleGlycan;

        typedef std::function<void(
            pcg32&,
            const AngleSettings& settings,
            WiggleGlycan,
            const assembly::Graph&,
            const AssemblyData&,
            MutableData&,
            const std::vector<GlycanShapePreference>&,
            const std::vector<size_t>&)>
            SidechainAdjustment;

        typedef std::function<GlycanShapePreference(
            pcg32& rng, const AngleSettings&, const AssemblyData&, const MutableData&, size_t glycanId)>
            GlycanShapeRandomizer;

        GlycoproteinState randomDescent(
            pcg32& rng,
            const AngleSettings& settings,
            WiggleGlycan wiggleGlycan,
            GlycanShapeRandomizer randomizeShape,
            SidechainAdjustment adjustSidechains,
            uint persistCycles,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            const GlycoproteinState& initialState);

        void resolveOverlapsWithWiggler(
            pcg32& rng,
            const AngleSettings& initialAngleSettings,
            const AngleSettings& mainAngleSettings,
            SidechainAdjustment adjustSidechains,
            SidechainAdjustment restoreSidechains,
            GlycanShapeRandomizer& randomizeShape,
            WiggleGlycan wiggleGlycan,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t persistCycles,
            bool deleteSitesUntilResolved);
    } // namespace gpbuilder
} // namespace gmml

#endif