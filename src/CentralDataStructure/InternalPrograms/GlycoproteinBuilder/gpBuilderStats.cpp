#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinUtil.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/structuredFiles.hpp"

#include <vector>
#include <string>
#include <sstream>

namespace glycoproteinBuilder
{
    Summary summarizeStats(const assembly::Graph& graph, const AssemblyData& data,
                           const GlycoproteinBuilderInputs& input, uint64_t seed,
                           const std::vector<StructureStats>& stats)
    {
        std::function<std::vector<std::string>(const StructureStats&)> toRow = [&](const StructureStats& stat)
        {
            auto glycanSitesToStream = [&input](std::ostream& stream, const std::vector<size_t>& glycans)
            {
                for (size_t n = 0; n < glycans.size(); n++)
                {
                    stream << input.glycositesInputVector[n].proteinResidueId << (n == glycans.size() - 1 ? "" : " ");
                }
            };
            std::ostringstream status;
            std::ostringstream nonviable;
            std::ostringstream deletions;
            std::ostringstream movedSidechainResidues;
            std::vector<size_t> movedSidechains = codeUtils::boolsToIndices(stat.residueSidechainMoved);
            if (stat.rejected)
            {
                status << "rejected";
            }
            else if (stat.deletions)
            {
                std::vector<size_t> deletedGlycans = codeUtils::boolsToIndices(codeUtils::indicesToValues(
                    codeUtils::vectorNot(stat.selection.molecules), data.glycans.moleculeId));
                status << "deletions";
                glycanSitesToStream(deletions, deletedGlycans);
            }
            std::vector<size_t> nonViableGlycans =
                codeUtils::boolsToIndices(codeUtils::indicesToValues(stat.nonViableMolecule, data.glycans.moleculeId));
            glycanSitesToStream(nonviable, nonViableGlycans);
            double highest = 0.0;
            for (size_t molecule : includedGlycanMoleculeIds(data, stat.selection.molecules))
            {
                for (size_t n : moleculeAtoms(graph, molecule))
                {
                    highest = std::max(highest, stat.overlap[n].weight);
                }
            }
            for (size_t n = 0; n < movedSidechains.size(); n++)
            {
                size_t residue = movedSidechains[n];
                movedSidechainResidues << data.residues.names[residue] << data.residues.numbers[residue]
                                       << (n == movedSidechains.size() - 1 ? "" : " ");
            }
            return std::vector<std::string> {stat.filename,   status.str(),    std::to_string(highest),
                                             nonviable.str(), deletions.str(), movedSidechainResidues.str()};
        };

        std::function<std::string(bool)> boolStr = [](bool b)
        {
            return (b ? "true" : "false");
        };

        codeUtils::TextTable parameterTable {
            {"Parameter", "Value"},
            {{numberOfSamplesParameter, std::to_string(input.numberOfSamples)},
             {seedParameter, std::to_string(seed)},
             {prepareForMDParameter, boolStr(input.MDprep)},
             {persistCyclesParameter, std::to_string(input.persistCycles)},
             {overlapRejectionThresholdParameter, std::to_string(input.overlapRejectionThreshold)},
             {useInitialGlycositeResidueConformationParameter, boolStr(input.useInitialGlycositeResidueConformation)},
             {moveOverlappingSidechainsParameter, boolStr(input.moveOverlappingSidechains)},
             {deleteIncompatibleSitesParameter, boolStr(input.deleteSitesUntilResolved)},
             {overlapToleranceParameter, std::to_string(input.overlapTolerance)},
             {ignoreHydrogenParameter, boolStr(input.ignoreHydrogen)}}
        };

        codeUtils::TextTable structureTable {
            {"Filename", "Status", "Highest Lennard-Jones potential", "Failed glycosites", "Deleted glycosites",
             "Moved sidechains"},
            codeUtils::vectorMap(toRow, stats)
        };

        return {input.inputFileName, input.substrateFileName, parameterTable, structureTable};
    }
} // namespace glycoproteinBuilder
