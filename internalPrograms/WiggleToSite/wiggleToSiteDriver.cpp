#include "include/internalPrograms/WiggleToSite/inputs.hpp"
#include "include/internalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/util/filesystem.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    using namespace gmml;
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile\n";
        std::cout << "Exmpl: " << argv[0] << " input.txt\n";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    std::cout << "Input file is " << inputFile << "\n";
    WiggleToSiteInputs inputStruct(inputFile);
    std::cout << "Reading input file complete\n" << std::flush;
    std::cout << inputStruct.Print();
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    const ParameterManager parameterManager = loadParameters(baseDir);
    WiggleToSite wiggler(parameterManager, inputStruct);
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
