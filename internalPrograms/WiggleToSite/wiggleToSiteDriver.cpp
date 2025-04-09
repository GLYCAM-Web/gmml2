#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/inputs.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile\n";
        std::cout << "Exmpl: " << argv[0] << " input.txt\n";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    std::cout << "Input file is " << inputFile << "\n";
    gmmlPrograms::WiggleToSiteInputs inputStruct(inputFile);
    std::cout << "Reading input file complete\n" << std::flush;
    std::cout << inputStruct.Print();
    std::string baseDir = codeUtils::toString(codeUtils::pathAboveCurrentExecutableDir());
    const cdsParameters::ParameterManager parameterManager(baseDir);
    gmmlPrograms::WiggleToSite wiggler(parameterManager, inputStruct);
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
