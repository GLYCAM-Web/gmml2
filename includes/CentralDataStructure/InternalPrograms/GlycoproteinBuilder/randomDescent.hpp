#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_RANDOMDESCENT_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapCount.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>

namespace glycoproteinBuilder
{
    struct GlycoproteinState
    {
        cds::Overlap overlap;
        std::vector<size_t> overlapSites;
        std::vector<std::vector<cds::ResidueLinkageShapePreference>> preferences;
    };

    typedef std::function<void(const AssemblyGraphs&, AssemblyData&, size_t glycanId, OverlapWeight,
                               std::vector<cds::ResidueLinkageShapePreference>&)>
        WiggleGlycan;

    typedef std::function<std::vector<cds::ResidueLinkageShapePreference>(const AssemblyGraphs&, const AssemblyData&,
                                                                          size_t glycanId)>
        GlycanShapeRandomizer;

    void wiggleLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId, size_t linkageId,
                       const cds::AngleSearchSettings& searchSettings, OverlapWeight weight,
                       const cds::ResidueLinkageShapePreference& shapePreference);
    void wiggleGlycan(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                      const cds::AngleSearchSettings& searchSettings, OverlapWeight weight,
                      const std::vector<cds::ResidueLinkageShapePreference>& preferences);
    GlycoproteinState randomDescent(pcg32 rng, GlycanShapeRandomizer randomizeShape,
                                    const cds::AngleSearchSettings& searchSettings, uint persistCycles,
                                    OverlapWeight overlapWeight, const AssemblyGraphs& graphs, AssemblyData& data,
                                    const GlycoproteinState& initialState);
} // namespace glycoproteinBuilder
#endif
