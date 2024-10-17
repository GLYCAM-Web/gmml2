#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>

namespace glycoproteinBuilder
{
    GlycoproteinBuilderInputs readGPInputFile(std::string inputFileName)
    {
        static const std::string proteinGlycanSection = "ProteinResidue, GlycanName";
        bool foundGlycanSection                       = false;
        bool readingGlycanSection                     = false;
        GlycoproteinBuilderInputs gpInputs;
        auto processLine = [&](const std::string& original, size_t lineNumber)
        {
            auto throwError = [&](const std::string& str)
            {
                throw std::runtime_error("Error reading input file on line " + std::to_string(lineNumber) + ": " + str +
                                         ":\n" + original + "\n");
            };
            auto parseUlong = [&](const std::string& str)
            {
                if (!std::all_of(str.begin(), str.end(), ::isdigit))
                {
                    throwError("'" + str + "' is not a valid non-negative integer");
                }
                try
                {
                    return std::stoul(str);
                }
                catch (...)
                {
                    throwError("'" + str + "' is not a valid non-negative integer");
                }
                return ulong(0);
            };
            auto parseBool = [&](const std::string& str)
            {
                if (str == "true")
                {
                    return true;
                }
                else if (str == "false")
                {
                    return false;
                }
                else
                {
                    throwError("'" + str + "' expected to be 'true' or 'false'");
                }
                return false;
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
                        gpInputs.glycositesInputVector.emplace_back(splitLine[0], splitLine[1]);
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
                    if (parameter == proteinGlycanSection)
                    {
                        foundGlycanSection   = true;
                        readingGlycanSection = true;
                        // return early to bypass checks below
                        return;
                    }
                    if (split.size() != 2)
                    {
                        throwError("input doesn't follow format 'parameter:value' or '" + proteinGlycanSection + ":'");
                    }
                    const std::string& value = split[1];
                    if (parameter == "Protein")
                    {
                        gpInputs.substrateFileName = value;
                    }
                    else if (parameter == "NumberOfOutputStructures")
                    {
                        gpInputs.number3DStructures = parseUlong(value);
                    }
                    else if (parameter == "persistCycles")
                    {
                        gpInputs.persistCycles = parseUlong(value);
                    }
                    else if (parameter == "freezeGlycositeResidueConformation")
                    {
                        gpInputs.freezeGlycositeResidueConformation = parseBool(value);
                    }
                    else if (parameter == "deleteIncompatibleSites")
                    {
                        gpInputs.deleteSitesUntilResolved = parseBool(value);
                    }
                    else if (parameter == "seed")
                    {
                        gpInputs.isDeterministic = true;
                        gpInputs.seed            = parseUlong(value);
                    }
                    else if (parameter == "skipMDPrep")
                    {
                        gpInputs.skipMDPrep = parseBool(value);
                    }
                    else if (parameter == "maxThreads")
                    {}
                    else if (parameter == "prepFileLocation")
                    {}
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
            throw std::runtime_error("Error reading input file: '" + proteinGlycanSection + ":' section missing\n");
        }
        if (readingGlycanSection)
        {
            throw std::runtime_error("Error reading input file: 'END' expected to close '" + proteinGlycanSection +
                                     ":' section'\n");
        }
        if (gpInputs.glycositesInputVector.empty())
        {
            throw std::runtime_error("Error reading input file: no glycosites requested\n");
        }
        return gpInputs;
    }
} // namespace glycoproteinBuilder
