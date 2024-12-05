#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/version.h"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"

#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,       "input-file", ""},
        {ArgReq::optional, ArgType::unnamed, "output-directory", ""},
        {ArgReq::optional,    ArgType::flag,           "--help", ""},
        {ArgReq::optional,    ArgType::flag,        "--version", ""}
    };
    codeUtils::Arguments arguments;
    try
    {
        arguments = codeUtils::readArguments(argc, argv);
        codeUtils::validateFlagsAndOptions(arguments, argumentDefinitions);
        std::string helpFlag = "--help";
        if (codeUtils::contains(arguments.names, helpFlag))
        {
            std::cout << codeUtils::helpString(arguments.programName, argumentDefinitions);
            std::cout << "\n"
                      << "For more information, see https://github.com/GLYCAM-Web/gmml2\n";
            std::exit(0);
        }
        std::string versionFlag = "--version";
        if (codeUtils::contains(arguments.names, versionFlag))
        {
            std::cout << "Glycoprotein Builder & GMML2 version " << GMML_VERSION << "\n";
            std::exit(0);
        }
        codeUtils::validateArgumentCount(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << codeUtils::helpString(arguments.programName, argumentDefinitions);
        std::exit(1);
    }

    try
    {
        std::string inputFile = arguments.unnamed[0];
        std::string outputDir = arguments.unnamed.size() > 1 ? arguments.unnamed[1] : ".";
        struct stat info;
        if (stat(outputDir.c_str(), &info) != 0)
        {
            std::cerr << "Folder " << outputDir << "/ does not exist and it isn't my job to make it.\n";
            std::exit(EXIT_FAILURE);
        }
        outputDir = outputDir + "/";
        std::cout << "Input file is " << inputFile << "\n";
        glycoproteinBuilder::GlycoproteinBuilderInputs inputStruct = glycoproteinBuilder::readGPInputFile(inputFile);
        std::cout << "Reading input file complete, on to construction\n" << std::flush;
        glycoproteinBuilder::GlycoproteinBuilder glycoproteinBuilder(inputStruct);
        std::cout << "Resolving overlaps" << std::endl;
        glycoproteinBuilder.ResolveOverlaps(outputDir);
    }
    catch (const std::runtime_error& error)
    {
        std::string errorMessage(error.what());
        std::string message = "glycoproteinBuilder: " + errorMessage;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        std::cout << message;
        std::exit(1);
    }
    catch (...)
    {
        std::string message =
            "glycoproteinBuilder: unexpected error. Please report how you got this to glycam@gmail.com";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
