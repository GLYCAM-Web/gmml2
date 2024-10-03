#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"

#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile outputDirectory\n";
        std::cout << "Exmpl: " << argv[0] << " input.txt output\n";
        std::exit(1);
    }
    try
    {
        std::string inputFile = argv[1];
        std::string outputDir = (argc >= 3) ? argv[2] : ".";
        struct stat info;
        if (stat(outputDir.c_str(), &info) != 0)
        {
            std::cerr << "Folder " << outputDir << "/ does not exist and it isn't my job to make it.\n";
            std::exit(EXIT_FAILURE);
        }
        outputDir = outputDir + "/";
        std::cout << "Test 017 Input file is " << inputFile << "\n";
        glycoproteinBuilder::GlycoproteinBuilderInputs inputStruct = glycoproteinBuilder::readGPInputFile(inputFile);
        std::cout << "Test 017 Reading input file complete, on to construction\n" << std::flush;
        glycoproteinBuilder::GlycoproteinBuilder glycoproteinBuilder(inputStruct);
        std::cout << "Test 017 Resolving overlaps" << std::endl;
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
    std::cout << "Test 017 Program got to end ok" << std::endl;
    return 0;
}
