#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/version.h"

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
    enum ARGUMENTS
    {
        SEQUENCE,
        FILENAME,
        HELP,
        VERSION,
        RELATIVE_PATHS
    };

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,       SEQUENCE,               "", ' ', "sequence"},
        {ArgReq::required, ArgType::unnamed,       FILENAME,               "", ' ', "filename"},
        {ArgReq::optional,    ArgType::flag,           HELP,           "help", 'h',         ""},
        {ArgReq::optional,    ArgType::flag,        VERSION,        "version", 'v',         ""},
        {ArgReq::optional,    ArgType::flag, RELATIVE_PATHS, "relative-paths", ' ',         ""}
    };
    std::string programName = codeUtils::programName(argv);
    codeUtils::Arguments arguments;
    try
    {
        arguments = codeUtils::readArguments(argc, argv, argumentDefinitions);
        if (codeUtils::contains<int>(arguments.ids, HELP))
        {
            std::cout << codeUtils::helpString(programName, argumentDefinitions);
            std::cout << "\n"
                      << "For more information, see https://github.com/GLYCAM-Web/gmml2\n";
            std::exit(0);
        }
        if (codeUtils::contains<int>(arguments.ids, VERSION))
        {
            std::cout << "GMML version " << GMML_VERSION << "\n";
            std::exit(0);
        }
        codeUtils::validateArguments(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << codeUtils::helpString(programName, argumentDefinitions);
        std::exit(1);
    }

    std::string baseDir  = codeUtils::toString(codeUtils::pathAboveCurrentExecutableDir());
    std::string SNFGDir  = codeUtils::SNFGSymbolsDir();
    std::string filename = "";
    std::string sequence;
    bool relative = false;

    for (const auto& arg : arguments.args)
    {
        switch (arg.id)
        {
            case ARGUMENTS::SEQUENCE:
                {
                    sequence = arg.value;
                    break;
                }
            case ARGUMENTS::FILENAME:
                {
                    filename = arg.value;
                    break;
                }
            case ARGUMENTS::RELATIVE_PATHS:
                {
                    relative = true;
                    break;
                }
            default:
                break;
        }
    }

    baseDir += "/";
    if (!relative)
    {
        SNFGDir = baseDir + SNFGDir;
        baseDir = "";
    }

    cdsCondensedSequence::GraphVizDotConfig config(baseDir, SNFGDir, filename);

    CondensedSequence::drawGlycan(config, sequence);
    return 0;
}
