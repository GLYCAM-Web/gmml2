#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/overlapCount.hpp"

#include <functional>

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
            const assembly::Graph& graph,
            const AssemblyData& data,
            MutableData& mutableData,
            size_t persistCycles,
            bool deleteSitesUntilResolved);
    } // namespace gpbuilder
} // namespace gmml

#endif