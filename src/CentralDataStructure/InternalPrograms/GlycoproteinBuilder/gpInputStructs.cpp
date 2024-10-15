#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>

using glycoprotein::GlycoproteinBuilderInputs;

GlycoproteinBuilderInputs glycoprotein::readGPInputFile(std::string inputFileName)
{
    bool readingGlycanSection = false;
    GlycoproteinBuilderInputs gpInputs;
    auto processLine = [&](const std::string& line)
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
                gpInputs.glycositesInputVector.emplace_back(splitLine.at(0), splitLine.at(1));
            }
        }
        else
        {
            if (codeUtils::startsWith(line, "ProteinResidue, GlycanName:"))
            {
                readingGlycanSection = true;
            }
            if (codeUtils::startsWith(line, "Protein:"))
            {
                gpInputs.substrateFileName = codeUtils::split(line, ':').at(1);
            }
            if (codeUtils::startsWith(line, "NumberOfOutputStructures:"))
            {
                gpInputs.number3DStructures = std::stoul(codeUtils::split(line, ':').at(1));
            }
            if (codeUtils::startsWith(line, "maxThreads:"))
            {
                gpInputs.maxThreads = std::stoul(codeUtils::split(line, ':').at(1));
            }
            if (codeUtils::startsWith(line, "persistCycles:"))
            {
                gpInputs.persistCycles = std::stoul(codeUtils::split(line, ':').at(1));
            }
            if (codeUtils::startsWith(line, "freezeGlycositeResidueConformation:"))
            {
                gpInputs.freezeGlycositeResidueConformation = codeUtils::split(line, ':').at(1) == "true";
            }
            if (codeUtils::startsWith(line, "seed:"))
            {
                gpInputs.isDeterministic = true;
                gpInputs.seed            = std::stoul(codeUtils::split(line, ':').at(1));
            }
            if (codeUtils::startsWith(line, "skipMDPrep:"))
            {
                gpInputs.skipMDPrep = codeUtils::split(line, ':').at(1) == "true";
            }
        }
    };
    codeUtils::readFileLineByLine(inputFileName, processLine);
    if (gpInputs.glycositesInputVector.empty())
    {
        throw std::runtime_error(
            "Error reading from gpInput file, no glycosites requested. Perhaps your formatting is incorrect.\n");
    }
    return gpInputs;
}
