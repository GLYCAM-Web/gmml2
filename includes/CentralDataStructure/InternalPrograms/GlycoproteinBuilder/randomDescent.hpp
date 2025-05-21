#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>

namespace glycoproteinBuilder
{
    struct GlycoproteinState
    {
        cds::Overlap overlap;
        std::vector<size_t> sitesWithOverlap;
        std::vector<size_t> sitesAboveOverlapThreshold;
        std::vector<cds::GlycanShapePreference> preferences;
    };

    typedef std::function<cds::GlycanShapePreference(pcg32& rng, const AngleSettings&, const AssemblyData&,
                                                     const MutableData&, size_t glycanId)>
        GlycanShapeRandomizer;

    GlycoproteinState randomDescent(pcg32& rng, const AngleSettings& settings, WiggleGlycan wiggleGlycan,
                                    GlycanShapeRandomizer randomizeShape, SidechainAdjustment adjustSidechains,
                                    uint persistCycles, const assembly::Graph& graph, const AssemblyData& data,
                                    MutableData& mutableData, const GlycoproteinState& initialState);
    void resolveOverlapsWithWiggler(pcg32& rng, const AngleSettings& initialAngleSettings,
                                    const AngleSettings& mainAngleSettings, SidechainAdjustment adjustSidechains,
                                    SidechainAdjustment restoreSidechains, GlycanShapeRandomizer& randomizeShape,
                                    WiggleGlycan wiggleGlycan, const assembly::Graph& graph, const AssemblyData& data,
                                    MutableData& mutableData, size_t persistCycles, bool deleteSitesUntilResolved);
} // namespace glycoproteinBuilder
#endif