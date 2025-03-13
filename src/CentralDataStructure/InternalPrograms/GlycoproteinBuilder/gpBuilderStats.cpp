#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/CodeUtils/containers.hpp"

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
            std::ostringstream status;
            std::ostringstream deletions;
            if (stat.rejected)
            {
                status << "rejected";
            }
            else if (stat.deletions)
            {
                std::vector<size_t> deletedGlycans = codeUtils::boolsToIndices(codeUtils::indicesToValues(
                    codeUtils::vectorNot(stat.selection.molecules), data.glycans.moleculeId));
                status << "deletions";
                for (size_t n = 0; n < deletedGlycans.size(); n++)
                {
                    deletions << input.glycositesInputVector[n].proteinResidueId
                              << (n == deletedGlycans.size() - 1 ? "" : " ");
                }
            }
            double highest = 0.0;
            for (size_t molecule : includedGlycanMoleculeIds(data, stat.selection.molecules))
            {
                for (size_t n : moleculeAtoms(graph, molecule))
                {
                    highest = std::max(highest, stat.overlap[n].weight);
                }
            }
            return std::vector<std::string> {stat.filename, status.str(), std::to_string(highest), deletions.str()};
        };

        std::function<std::string(bool)> boolStr = [](bool b)
        {
            return (b ? "true" : "false");
        };

        Table parameterTable {
            {"Parameter", "Value"},
            {{numberOfStructuresParameter, std::to_string(input.number3DStructures)},
             {seedParameter, std::to_string(seed)},
             {prepareForMDParameter, boolStr(input.MDprep)},
             {persistCyclesParameter, std::to_string(input.persistCycles)},
             {overlapRejectionThresholdParameter,
              (input.rejectExcessiveGlycanOverlaps ? std::to_string(input.overlapRejectionThreshold) : "none")},
             {useInitialGlycositeResidueConformationParameter, boolStr(input.useInitialGlycositeResidueConformation)},
             {moveOverlappingSidechainsParameter, boolStr(input.moveOverlappingSidechains)},
             {deleteIncompatibleSitesParameter, boolStr(input.deleteSitesUntilResolved)}}
        };

        Table structureTable {
            {"Filename", "Status", "Highest atom potential", "Deleted sites"},
            codeUtils::vectorMap(toRow, stats)
        };

        return {input.inputFileName, input.substrateFileName, parameterTable, structureTable};
    }

    std::string htmlSummary(const Summary& summary, const std::vector<std::string>& headerLines)
    {
        auto printTable = [](const Table& table)
        {
            std::ostringstream ss;
            auto printRow = [&](const std::vector<std::string>& row)
            {
                ss << "<tr>\n";
                for (auto& str : row)
                {
                    ss << "<th>" << str << "</th>\n";
                }
                ss << "</tr>\n";
            };
            ss << "<table>\n";
            printRow(table.header);
            for (auto& row : table.rows)
            {
                printRow(row);
            }
            ss << "</table>\n";
            return ss.str();
        };
        std::ostringstream ss;

        ss << "<!DOCTYPE html>\n"
           << "<html>\n"
           << "<head>\n"
           << "<style>\n"
           << "table {\n"
           << "font-family: arial, sans-serif;\n"
           << "border-collapse: collapse;\n"
           << "}\n"
           << "td, th {\n"
           << "border: 1px solid #dddddd;\n"
           << "text-align: left;\n"
           << "padding: 8px;\n"
           << "}\n"
           << "tr:nth-child(even) {\n"
           << "background-color: #dddddd;\n"
           << "}\n"
           << "</style>\n"
           << "</head>\n"
           << "<body>\n";

        ss << "<h2>Glycoprotein Builder</h2>\n";

        auto printParagraph = [](std::ostringstream& ss, const std::string& str)
        {
            ss << "<p>" << str << "</p>\n";
        };

        for (auto& line : headerLines)
        {
            printParagraph(ss, line);
        }

        ss << "<h3>Input</h3>\n";
        ss << "<p>Filename: " << summary.filename << "</p>\n";
        ss << "<p>Protein: " << summary.proteinFilename << "</p>\n";
        ss << printTable(summary.parameterTable);

        ss << "<h3>Structures</h3>\n";
        ss << printTable(summary.structuretable);

        ss << "</body>\n"
           << "</html>\n";

        return ss.str();
    }

    std::string plaintextSummary(const Summary& summary, const std::vector<std::string>& headerLines)
    {
        auto printTable = [](const Table& table)
        {
            size_t columns = table.header.size();
            std::vector<size_t> columnSize(columns, 0);
            for (size_t n = 0; n < columns; n++)
            {
                size_t size = table.header[n].size();
                for (auto& row : table.rows)
                {
                    size = std::max(size, row[n].size());
                }
                columnSize[n] = size;
            }
            std::ostringstream ss;
            auto printRow = [&](const std::string& delimiter, const std::vector<std::string>& row)
            {
                for (size_t n = 0; n < columns; n++)
                {
                    const std::string& str = row[n];
                    ss << str << std::string(columnSize[n] - str.size(), ' ');
                    if (n < columns - 1)
                    {
                        ss << delimiter;
                    }
                }
                ss << "\n";
            };
            printRow(" | ", table.header);
            for (auto& row : table.rows)
            {
                printRow("   ", row);
            }
            return ss.str();
        };
        std::ostringstream ss;

        ss << "### Glycoprotein Builder ###\n";

        auto printParagraph = [](std::ostringstream& ss, const std::string& str)
        {
            ss << str << "\n";
        };

        for (auto& line : headerLines)
        {
            printParagraph(ss, line);
        }

        ss << "\n";
        ss << "## Input ##\n";
        ss << "Filename: " << summary.filename << "\n";
        ss << "Protein: " << summary.proteinFilename << "\n";
        ss << "\n";
        ss << printTable(summary.parameterTable);
        ss << "\n";

        ss << "## Structures ##\n";
        ss << "\n";
        ss << printTable(summary.structuretable);

        return ss.str();
    }

    std::string csvTable(const Table& table)
    {
        std::ostringstream ss;
        auto printRow = [&](const std::string& linebreak, const std::vector<std::string>& row)
        {
            for (size_t n = 0; n < row.size(); n++)
            {
                ss << row[n] << (n < row.size() - 1 ? "," : linebreak);
            }
        };
        printRow("\n", table.header);
        for (size_t n = 0; n < table.rows.size(); n++)
        {
            std::string linebreak = (n == table.rows.size() - 1) ? "" : "\n";
            printRow(linebreak, table.rows[n]);
        }
        return ss.str();
    }
} // namespace glycoproteinBuilder
