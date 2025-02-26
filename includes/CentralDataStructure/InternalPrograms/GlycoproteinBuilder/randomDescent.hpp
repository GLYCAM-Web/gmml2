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
        std::vector<size_t> overlapSites;
        std::vector<cds::GlycanShapePreference> preferences;
    };

    typedef std::function<void(const assembly::Graph&, const AssemblyData&, size_t glycanId, const OverlapMultiplier&,
                               cds::GlycanShapePreference&)>
        WiggleGlycan;

    typedef std::function<cds::GlycanShapePreference(pcg32& rng, const AssemblyData&, const MutableData&,
                                                     size_t glycanId)>
        GlycanShapeRandomizer;

    void wiggleLinkage(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                       const std::vector<bool>& includedAtoms, size_t glycanId, size_t linkageId,
                       const cds::AngleSearchSettings& searchSettings, const OverlapMultiplier& overlapMultiplier,
                       const cds::ResidueLinkageShapePreference& shapePreference);
    void wiggleGlycan(const assembly::Graph& graph, const AssemblyData& data, MutableData& mutableData,
                      const std::vector<bool>& includedAtoms, size_t glycanId,
                      const cds::AngleSearchSettings& searchSettings, const OverlapMultiplier& overlapMultiplier,
                      const cds::GlycanShapePreference& preferences);
    GlycoproteinState randomDescent(pcg32& rng, GlycanShapeRandomizer randomizeShape,
                                    SidechainAdjustment adjustSidechains,
                                    const cds::AngleSearchSettings& searchSettings, uint persistCycles,
                                    const OverlapMultiplier& overlapMultiplier, const assembly::Graph& graph,
                                    const AssemblyData& data, MutableData& mutableData,
                                    const GlycoproteinState& initialState);
} // namespace glycoproteinBuilder
#endif