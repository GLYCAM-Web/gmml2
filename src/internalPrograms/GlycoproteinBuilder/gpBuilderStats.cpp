#include "include/internalPrograms/GlycoproteinBuilder/gpBuilderStats.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinUtil.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/util/containers.hpp"
#include "include/util/structuredFiles.hpp"

#include <sstream>
#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        Summary summarizeStats(
            const assembly::Graph& graph,
            const AssemblyData& data,
            const GlycoproteinBuilderInputs& input,
            uint64_t seed,
            const std::vector<StructureStats>& stats)
        {
            std::function<std::vector<std::string>(const StructureStats&)> toRow = [&](const StructureStats& stat)
            {
                auto glycanSitesToStream = [&input](std::ostream& stream, const std::vector<size_t>& glycans)
                {
                    for (size_t n = 0; n < glycans.size(); n++)
                    {
                        stream << input.glycositesInputVector[n].proteinResidueId
                               << (n == glycans.size() - 1 ? "" : " ");
                    }
                };
                std::ostringstream status;
                std::ostringstream nonviable;
                std::ostringstream deletions;
                std::ostringstream movedSidechainResidues;
                std::vector<size_t> movedSidechains = util::boolsToIndices(stat.residueSidechainMoved);
                if (stat.rejected)
                {
                    status << "rejected";
                }
                else if (stat.deletions)
                {
                    std::vector<size_t> deletedGlycans = util::boolsToIndices(
                        util::indicesToValues(util::vectorNot(stat.selection.molecules), data.glycans.moleculeId));
                    status << "deletions";
                    glycanSitesToStream(deletions, deletedGlycans);
                }
                std::vector<size_t> nonViableGlycans =
                    util::boolsToIndices(util::indicesToValues(stat.nonViableMolecule, data.glycans.moleculeId));
                glycanSitesToStream(nonviable, nonViableGlycans);
                double highest = 0.0;
                for (size_t molecule : includedGlycanMoleculeIds(data, stat.selection.molecules))
                {
                    for (size_t n : moleculeAtoms(graph, molecule))
                    {
                        highest = std::max(highest, stat.overlap[n]);
                    }
                }
                for (size_t n = 0; n < movedSidechains.size(); n++)
                {
                    size_t residue = movedSidechains[n];
                    movedSidechainResidues << data.residues.names[residue] << data.residues.numbers[residue]
                                           << (n == movedSidechains.size() - 1 ? "" : " ");
                }
                return std::vector<std::string> {
                    stat.filename,
                    status.str(),
                    std::to_string(highest),
                    nonviable.str(),
                    deletions.str(),
                    movedSidechainResidues.str()};
            };

            std::function<std::string(bool)> boolStr = [](bool b) { return (b ? "true" : "false"); };

            util::TextTable parameterTable {
                {"Parameter", "Value"},
                {{numberOfSamplesParameter, std::to_string(input.numberOfSamples)},
                 {seedParameter, std::to_string(seed)},
                 {prepareForMDParameter, boolStr(input.MDprep)},
                 {persistCyclesParameter, std::to_string(input.persistCycles)},
                 {overlapRejectionThresholdParameter, std::to_string(input.overlapRejectionThreshold)},
                 {useInitialGlycositeResidueConformationParameter,
                  boolStr(input.useInitialGlycositeResidueConformation)},
                 {moveOverlappingSidechainsParameter, boolStr(input.moveOverlappingSidechains)},
                 {deleteIncompatibleSitesParameter, boolStr(input.deleteSitesUntilResolved)},
                 {overlapToleranceParameter, std::to_string(input.overlapTolerance)},
                 {ignoreHydrogenParameter, boolStr(input.ignoreHydrogen)}}
            };

            util::TextTable structureTable {
                {"Filename",
                 "Status", "Highest Lennard-Jones potential",
                 "Failed glycosites", "Deleted glycosites",
                 "Moved sidechains"},
                util::vectorMap(toRow, stats)
            };

            return {input.inputFileName, input.substrateFileName, parameterTable, structureTable};
        }
    } // namespace gpbuilder
} // namespace gmml
