#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANWIGGLE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCANWIGGLE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

namespace glycoproteinBuilder
{
    void wiggleLinkage(const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection,
                       MutableData& mutableData, size_t linkageId, const cds::AngleSearchSettings& searchSettings,
                       const cds::ResidueLinkageShapePreference& shapePreference);
    void wiggleGlycan(const assembly::Graph& graph, const AssemblyData& data, const assembly::Selection& selection,
                      const cds::AngleSearchSettings& searchSettings, const cds::GlycanShapePreference& preferences,
                      MutableData& mutableData, size_t glycanId);
} // namespace glycoproteinBuilder
#endif
