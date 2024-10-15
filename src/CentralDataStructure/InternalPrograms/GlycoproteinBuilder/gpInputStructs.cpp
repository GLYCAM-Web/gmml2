#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <algorithm> // remove

using glycoprotein::GlycoproteinBuilderInputs;

GlycoproteinBuilderInputs glycoprotein::readGPInputFile(std::string inputFileName)
{
    //    std::cout << "About to read " << inputFileName << std::endl << std::flush;
    std::ifstream infile(inputFileName);
    if (!infile)
    {
        std::string message = "Uh oh, input file: " + inputFileName + ", could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    GlycoproteinBuilderInputs gpInputs;
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        strInput.erase(std::remove(strInput.begin(), strInput.end(), '\r'),
                       strInput.end()); // files created in windows, being read in unix.
        gmml::log(__LINE__, __FILE__, gmml::INF, strInput);
        if (codeUtils::startsWith(strInput, "Protein:"))
        {
            gpInputs.substrateFileName = codeUtils::split(strInput, ':').at(1);
        }
        if (codeUtils::startsWith(strInput, "NumberOfOutputStructures:"))
        {
            gpInputs.number3DStructures = std::stoul(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "maxThreads:"))
        {
            gpInputs.maxThreads = std::stoul(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "persistCycles:"))
        {
            gpInputs.persistCycles = std::stoul(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "freezeGlycositeResidueConformation:"))
        {
            gpInputs.freezeGlycositeResidueConformation = codeUtils::split(strInput, ':').at(1) == "true";
        }
        if (codeUtils::startsWith(strInput, "seed:"))
        {
            gpInputs.isDeterministic = true;
            gpInputs.seed            = std::stoul(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "skipMDPrep:"))
        {
            gpInputs.skipMDPrep = codeUtils::split(strInput, ':').at(1) == "true";
        }
        if (codeUtils::startsWith(strInput, "ProteinResidue, GlycanName:"))
        {
            std::string tempBuffer; //  Temporarily holds whatever getline() finds on the line;
            while (getline(infile, tempBuffer))
            {
                tempBuffer.erase(std::remove(tempBuffer.begin(), tempBuffer.end(), '\r'),
                                 tempBuffer.end()); // files created in windows, being read in unix.
                if (tempBuffer == "END")
                {
                    break;
                }
                std::vector<std::string> splitLine = codeUtils::split(tempBuffer, '|');
                gpInputs.glycositesInputVector.emplace_back(splitLine.at(0), splitLine.at(1));
            }
        }
    }
    //    std::cout << "Reading input file complete, just making a quick check\n" << std::flush;
    if (gpInputs.glycositesInputVector.empty())
    {
        throw std::runtime_error(
            "Error reading from gpInput file, no glycosites requested. Perhaps your formatting is incorrect.\n");
    }
    return gpInputs;
}
