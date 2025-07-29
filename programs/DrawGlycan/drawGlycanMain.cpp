#include "include/sequence/graphViz.hpp"
#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/sequence/sequencePrinter.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containers.hpp"
#include "include/util/filesystem.hpp"
#include "include/version.h"

#include <iostream>
#include <string>
#include <vector>

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

    using namespace gmml;
    using namespace gmml::sequence;
    using util::ArgReq;
    using util::ArgType;
    std::vector<util::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,       SEQUENCE,               "", ' ', "sequence"},
        {ArgReq::required, ArgType::unnamed,       FILENAME,               "", ' ', "filename"},
        {ArgReq::optional,    ArgType::flag,           HELP,           "help", 'h',         ""},
        {ArgReq::optional,    ArgType::flag,        VERSION,        "version", 'v',         ""},
        {ArgReq::optional,    ArgType::flag, RELATIVE_PATHS, "relative-paths", ' ',         ""}
    };
    std::string programName = util::programName(argv);
    util::Arguments arguments;
    try
    {
        arguments = util::readArguments(argc, argv, argumentDefinitions);
        if (util::contains<int>(arguments.ids, HELP))
        {
            std::cout << util::helpString(programName, argumentDefinitions);
            std::cout << "\n"
                      << "For more information, see https://github.com/GLYCAM-Web/gmml2\n";
            std::exit(0);
        }
        if (util::contains<int>(arguments.ids, VERSION))
        {
            std::cout << "GMML2 version " << GMML_VERSION << "\n";
            std::exit(0);
        }
        util::validateArguments(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << util::helpString(programName, argumentDefinitions);
        std::exit(1);
    }

    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    std::string SNFGDir = util::SNFGSymbolsDir();
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

    sequence::GraphVizDotConfig config(baseDir, SNFGDir, filename);

    printGraphViz(config, instantiate(parseSequence(sequence)));
    return 0;
}
