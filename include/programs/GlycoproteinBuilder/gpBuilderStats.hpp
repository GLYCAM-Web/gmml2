#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GPBUILDERSTATS_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GPBUILDERSTATS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/programs/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/programs/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/util/structuredFiles.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct StructureStats
        {
            std::string filename;
            bool rejected;
            bool deletions;
            assembly::Bounds bounds;
            assembly::Selection selection;
            std::vector<bool> nonViableMolecule;
            std::vector<bool> residueSidechainMoved;
            std::vector<double> overlap;
        };

        struct Summary
        {
            std::string filename;
            std::string proteinFilename;
            util::TextTable parameterTable;
            util::TextTable structuretable;
        };

        Summary summarizeStats(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const GlycoproteinBuilderInputs& input,
            uint64_t seed,
            const std::vector<StructureStats>& stats);
    } // namespace gpbuilder
} // namespace gmml

#endif
