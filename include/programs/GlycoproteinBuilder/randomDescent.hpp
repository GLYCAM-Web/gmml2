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
            double totalOverlap;
            OverlapSites overlapSites;
            std::vector<GlycanShapePreference> preferences;
            MutableData mutableData;
        };

        typedef std::function<AngleSettings(uint)> PersistCycleAngleSettings;

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
            pcg32&, const assembly::Graph&, const AssemblyData&, MutableData&, const std::vector<size_t>&)>
            SidechainAdjustment;

        typedef std::function<GlycanShapePreference(
            pcg32& rng,
            const AngleSettings&,
            const AssemblyData&,
            const assembly::Bounds&,
            const LinkageShapeSettings&,
            const std::vector<GlycanShapePreference>&,
            size_t)>
            GlycanShapeRandomizer;

        typedef std::function<GlycanShapePreference(
            pcg32& rng, const AssemblyData&, const assembly::Bounds&, size_t glycanId)>
            FullGlycanDihedralRandomizer;

        GlycoproteinState randomDescent(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            PersistCycleAngleSettings toAngleSettings,
            GlycanShapeRandomizer randomizeShape,
            const LinkageShapeSettings& shapeSettings,
            SidechainAdjustment adjustSidechains,
            uint persistCycles,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const GlycoproteinState& initialState);

        std::vector<GlycanShapePreference> randomizeInitialShapePreference(
            pcg32& rng,
            const AngleSettings& angleSettings,
            GlycanShapeRandomizer& randomizeShape,
            const LinkageShapeSettings& shapeSettings,
            const AssemblyData& data,
            const assembly::Bounds& bounds);

        GlycoproteinState resolveOverlapsWithWiggler(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            PersistCycleAngleSettings toAngleSettings,
            SidechainAdjustment adjustSidechains,
            SidechainAdjustment restoreSidechains,
            GlycanShapeRandomizer& randomizeShape,
            FullGlycanDihedralRandomizer& randomizeAllDihedrals,
            const LinkageShapeSettings& shapeSettings,
            const std::vector<GlycanShapePreference>& initialPreference,
            const OverlapSettings& overlapSettings,
            const assembly::Graph& graph,
            const AssemblyData& data,
            const MutableData& initialState,
            size_t persistCycles,
            size_t randomizationCycles,
            bool deleteSitesUntilResolved);
    } // namespace gpbuilder
} // namespace gmml

#endif