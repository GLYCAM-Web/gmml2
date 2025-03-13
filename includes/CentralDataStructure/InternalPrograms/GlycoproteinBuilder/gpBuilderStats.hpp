#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPBUILDERSTATS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPBUILDERSTATS_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Assembly/assemblyBounds.hpp"

#include <vector>
#include <string>

namespace glycoproteinBuilder
{
    struct StructureStats
    {
        std::string filename;
        bool rejected;
        bool deletions;
        assembly::Bounds bounds;
        assembly::Selection selection;
        std::vector<bool> residueSidechainMoved;
        std::vector<cds::Overlap> overlap;
    };

    struct Table
    {
        std::vector<std::string> header;
        std::vector<std::vector<std::string>> rows;
    };

    struct Summary
    {
        std::string filename;
        std::string proteinFilename;
        Table parameterTable;
        Table structuretable;
    };

    Summary summarizeStats(const assembly::Graph& graph, const AssemblyData& data,
                           const GlycoproteinBuilderInputs& input, uint64_t seed,
                           const std::vector<StructureStats>& stats);
    std::string htmlSummary(const Summary& summary, const std::vector<std::string>& headerLines);
    std::string plaintextSummary(const Summary& summary, const std::vector<std::string>& headerLines);
    std::string csvTable(const Table& table);
} // namespace glycoproteinBuilder
#endif
