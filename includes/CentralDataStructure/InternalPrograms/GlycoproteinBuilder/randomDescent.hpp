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

    void wiggleLinkage(const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection,
                       MutableData& mutableData, size_t linkageId, const cds::AngleSearchSettings& searchSettings,
                       const cds::ResidueLinkageShapePreference& shapePreference);
    void wiggleGlycan(const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection,
                      const cds::AngleSearchSettings& searchSettings, const cds::GlycanShapePreference& preferences,
                      MutableData& mutableData, size_t glycanId);
    GlycoproteinState randomDescent(pcg32& rng, const AngleSettings& settings, WiggleGlycan wiggleGlycan,
                                    GlycanShapeRandomizer randomizeShape, SidechainAdjustment adjustSidechains,
                                    uint persistCycles, const assembly::Graph& graph, const AssemblyData& data,
                                    MutableData& mutableData, const GlycoproteinState& initialState);
} // namespace glycoproteinBuilder
#endif