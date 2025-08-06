#include "include/carbohydrate/carbohydrate.hpp"
#include "include/carbohydrate/parameterManager.hpp"
#include "include/fileType/dot/graphViz.hpp"
#include "include/sequence/sequenceManipulation.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/sequence/sequencePrinter.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"
#include "include/version.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    enum ARGUMENTS
    {
        INPUT_FILE,
        DELIMITER,
        OUTPUT_DIR,
        HELP,
        VERSION,
        TEST_MODE,
        OVERWRITE_EXISTING
    };

    using namespace gmml;
    using namespace sequence;
    using util::ArgReq;
    using util::ArgType;
    const std::string overwriteFlag = "overwrite-existing-files";
    std::vector<util::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,         INPUT_FILE,            "", ' ',       "input-file"},
        {ArgReq::required, ArgType::unnamed,          DELIMITER,            "", ' ',   "list-delimiter"},
        {ArgReq::required, ArgType::unnamed,         OUTPUT_DIR,            "", ' ', "output-directory"},
        {ArgReq::optional,    ArgType::flag,               HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,            VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,          TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,    ArgType::flag, OVERWRITE_EXISTING, overwriteFlag, ' ',                 ""}
    };

    std::string programName = util::programName(argv);
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    std::string SNFGDir = util::SNFGSymbolsDir();

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
            std::cout << "Carbohydrate Builder & GMML2 version " << GMML_VERSION << "\n";
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

    try
    {
        std::string inputFile = "";
        char delimiter = '_';
        std::string outputDir = ".";
        std::string headerBaseString = "Produced by GMML2 (https://github.com/GLYCAM-Web/gmml2)";
        std::vector<std::string> headerLines {headerBaseString + " version " + std::string(GMML_VERSION)};
        bool testMode = false;
        bool overwrite = false;
        for (const auto& arg : arguments.args)
        {
            switch (arg.id)
            {
                case ARGUMENTS::INPUT_FILE:
                    {
                        inputFile = arg.value;
                        break;
                    }
                case ARGUMENTS::DELIMITER:
                    {
                        if (arg.value.size() != 1)
                        {
                            std::cout << "list delimiter must be a single character\n";
                            std::cout << util::helpString(programName, argumentDefinitions);
                            std::exit(1);
                        }
                        delimiter = arg.value[0];
                        break;
                    }
                case ARGUMENTS::OUTPUT_DIR:
                    {
                        outputDir = arg.value;
                        util::createDirectories(outputDir);
                        break;
                    }
                case ARGUMENTS::TEST_MODE:
                    {
                        testMode = true;
                        headerLines = {headerBaseString + " in test mode"};
                        std::cout << "Running in test mode\n";
                        break;
                    }
                case ARGUMENTS::OVERWRITE_EXISTING:
                    {
                        overwrite = true;
                        break;
                    }
                default:
                    break;
            }
        }
        if (!(testMode || overwrite || outputDir == "."))
        {
            if (!util::directoryIsEmptyOrNonexistent(outputDir))
            {
                throw std::runtime_error(
                    "Output directory '" + outputDir +
                    "' not empty. Please empty/remove it or run the program with the flag --" + overwriteFlag);
            }
        }

        std::string dotBaseDir = baseDir + "/";
        if (!testMode)
        {
            SNFGDir = dotBaseDir + SNFGDir;
            dotBaseDir = "";
        }

        struct SequenceInput
        {
            std::string id;
            std::string sequence;
        };

        std::vector<SequenceInput> lines;
        auto processLine = [&](const std::string& line, size_t)
        {
            if (!line.empty())
            {
                std::vector<std::string> splitLine = util::split(line, delimiter);
                if (splitLine.size() == 2)
                {
                    lines.push_back({splitLine[0], splitLine[1]});
                }
                else
                {
                    std::cerr << "Encountered problem when splitting this line >>>" << line << "<<< from file >>>"
                              << argv[1] << "<<< into an ID and carb string separated by your specified delimiter: >>>"
                              << argv[2] << "<<<\n";
                    std::exit(EXIT_FAILURE);
                }
            }
        };
        util::readFileLineByLine(inputFile, processLine);
        const ParameterManager parameterManager = loadParameters(baseDir);
        const util::SparseVector<double> elementRadii = vanDerWaalsRadii();
        for (auto& line : lines)
        {
            try
            {
                std::cout << "\n*********************\nBuilding " << line.id << ": " << line.sequence
                          << "\n*********************\n";
                SequenceData sequenceData = parseAndReorder(line.sequence);
                Molecule carbohydrate;
                std::vector<ResidueLinkage> linkages;
                initializeCarbohydrate(carbohydrate, linkages, parameterManager, elementRadii, sequenceData);
                generate3DStructureFiles(carbohydrate, outputDir, line.id, headerLines);
                dot::Config config {dotBaseDir, SNFGDir, outputDir + "/" + line.id + ".dot"};
                printGraphViz(config, instantiate(parseSequence(line.sequence)));
            }
            catch (const std::runtime_error& error)
            {
                util::log(__LINE__, __FILE__, util::ERR, error.what());
                std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << error.what()
                          << std::endl;
            }
            catch (...)
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::ERR,
                    "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
                std::cerr
                    << "ERROR carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report "
                       "how you got to this to glycam@gmail.com.";
            }
        }
    }
    catch (const std::runtime_error& error)
    {
        std::string errorMessage(error.what());
        std::string message = "carbohydrateBuilder: " + errorMessage;
        util::log(__LINE__, __FILE__, util::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
    catch (...)
    {
        std::string message =
            "carbohydrateBuilder: unexpected error. Please report how you got this to glycam@gmail.com";
        util::log(__LINE__, __FILE__, util::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
}
