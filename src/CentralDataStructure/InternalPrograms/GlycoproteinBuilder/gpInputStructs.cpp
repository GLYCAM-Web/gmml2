#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/parsing.hpp"

#include <vector>
#include <string>
#include <optional>
#include <stdexcept>

namespace glycoproteinBuilder
{
    GlycoproteinBuilderInputs readGPInputFile(std::string inputFileName)
    {
        const std::string glycanSectionParameter          = "ProteinResidue, GlycanName";
        const std::vector<std::string> requiredParameters = {proteinParameter};
        std::vector<std::string> foundParameters          = {};

        bool foundGlycanSection   = false;
        bool readingGlycanSection = false;
        GlycoproteinBuilderInputs gpInputs;
        gpInputs.inputFileName = inputFileName;
        auto processLine       = [&](const std::string& original, size_t lineNumber)
        {
            auto throwError = [&](const std::string& str)
            {
                throw std::runtime_error("Error reading input file on line " + std::to_string(lineNumber) + ": " + str +
                                         ":\n" + original + "\n");
            };
            auto parseUlong = [&](const std::string& str)
            {
                std::optional<ulong> opt = codeUtils::parseUlong(str);
                if (!opt.has_value())
                {
                    throwError("'" + str + "' is not a valid non-negative integer");
                }
                return opt.value();
            };
            auto parseDouble = [&](const std::string& str)
            {
                std::optional<double> opt = codeUtils::parseDouble(str);
                if (!opt.has_value())
                {
                    throwError("'" + str + "' is not a valid floating point number");
                }
                return opt.value();
            };
            auto parseBool = [&](const std::string& str)
            {
                std::optional<bool> opt = codeUtils::parseBool(str);
                if (!opt.has_value())
                {
                    throwError("'" + str + "' expected to be 'true' or 'false'");
                }
                return opt.value();
            };
            const std::string line = codeUtils::trimmedOfWhitespace(original);
            if (!line.empty())
            {
                if (readingGlycanSection)
                {
                    if (line == "END")
                    {
                        readingGlycanSection = false;
                    }
                    else
                    {
                        std::vector<std::string> splitLine = codeUtils::split(line, '|');
                        if (splitLine.size() != 2)
                        {
                            throwError("input doesn't follow format 'ProteinResidue|GlycanName' or 'END'");
                        }
                        gpInputs.glycositesInputVector.push_back({splitLine[0], splitLine[1]});
                    }
                }
                else
                {
                    const std::vector<std::string> split = codeUtils::split(line, ':');
                    if (line.find(":") == std::string::npos)
                    {
                        throwError("input doesn't follow format 'parameter:value'");
                    }
                    const std::string& parameter = split[0];
                    if (codeUtils::contains(foundParameters, parameter))
                    {
                        throwError("duplicate parameter '" + parameter + "'");
                    }
                    foundParameters.push_back(parameter);
                    if (parameter == glycanSectionParameter)
                    {
                        foundGlycanSection   = true;
                        readingGlycanSection = true;
                        // return early to bypass checks below
                        return;
                    }
                    if (split.size() != 2)
                    {
                        throwError("input doesn't follow format 'parameter:value' or '" + glycanSectionParameter +
                                   ":'");
                    }
                    const std::string& value = split[1];
                    if (parameter == proteinParameter)
                    {
                        gpInputs.substrateFileName = value;
                    }
                    else if (parameter == numberOfStructuresParameter)
                    {
                        gpInputs.number3DStructures = parseUlong(value);
                    }
                    else if (parameter == persistCyclesParameter)
                    {
                        gpInputs.persistCycles = parseUlong(value);
                    }
                    else if (parameter == overlapToleranceParameter)
                    {
                        gpInputs.overlapTolerance = parseDouble(value);
                    }
                    else if (parameter == freezeGlycositeResidueConformationParameter)
                    {
                        gpInputs.freezeGlycositeResidueConformation = parseBool(value);
                    }
                    else if (parameter == allowSidechainAdjustmentParameter)
                    {
                        gpInputs.allowSidechainAdjustment = parseBool(value);
                    }
                    else if (parameter == deleteIncompatibleSitesParameter)
                    {
                        gpInputs.deleteSitesUntilResolved = parseBool(value);
                    }
                    else if (parameter == atomPotentialRejectionThresholdParameter)
                    {
                        gpInputs.rejectExcessiveGlycanOverlaps   = true;
                        gpInputs.atomPotentialRejectionThreshold = parseDouble(value);
                    }
                    else if (parameter == seedParameter)
                    {
                        gpInputs.isDeterministic = true;
                        gpInputs.seed            = parseUlong(value);
                    }
                    else if (parameter == prepareForMDParameter)
                    {
                        bool b                = parseBool(value);
                        gpInputs.MDprep       = b;
                        gpInputs.writeOffFile = b;
                    }
                    else
                    {
                        throwError("unknown parameter '" + parameter + "'");
                    }
                }
            }
        };
        codeUtils::readFileLineByLine(inputFileName, processLine);
        if (!foundGlycanSection)
        {
            throw std::runtime_error("Error reading input file: '" + glycanSectionParameter + ":' section missing\n");
        }
        if (readingGlycanSection)
        {
            throw std::runtime_error("Error reading input file: 'END' expected to close '" + glycanSectionParameter +
                                     ":' section'\n");
        }
        for (auto& req : requiredParameters)
        {
            if (!codeUtils::contains(foundParameters, req))
            {
                throw std::runtime_error("Error reading input file: required parameter '" + req + "' not found\n");
            }
        }
        if (gpInputs.glycositesInputVector.empty())
        {
            throw std::runtime_error("Error reading input file: no glycosites requested\n");
        }
        return gpInputs;
    }
} // namespace glycoproteinBuilder
