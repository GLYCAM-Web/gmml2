#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/version.h"

#include <filesystem>
#include <iostream>
#include <string>
#include <fstream>

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

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    const std::string overwriteFlag                    = "overwrite-existing-files";
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,         INPUT_FILE,            "", ' ',       "input-file"},
        {ArgReq::required, ArgType::unnamed,          DELIMITER,            "", ' ',   "list-delimiter"},
        {ArgReq::required, ArgType::unnamed,         OUTPUT_DIR,            "", ' ', "output-directory"},
        {ArgReq::optional,    ArgType::flag,               HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,            VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,          TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,    ArgType::flag, OVERWRITE_EXISTING, overwriteFlag, ' ',                 ""}
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
            std::cout << "Carbohydrate Builder & GMML version " << GMML_VERSION << "\n";
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

    try
    {
        std::string inputFile        = "";
        char delimiter               = '_';
        std::string outputDir        = ".";
        std::string headerBaseString = "Produced by GMML (https://github.com/GLYCAM-Web/gmml2)";
        std::vector<std::string> headerLines {headerBaseString + " version " + std::string(GMML_VERSION)};
        int numThreads = 1;
        bool testMode  = false;
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
                            std::cout << codeUtils::helpString(programName, argumentDefinitions);
                            std::exit(1);
                        }
                        delimiter = arg.value[0];
                    }
                case ARGUMENTS::OUTPUT_DIR:
                    {
                        outputDir = arg.value;
                        codeUtils::createDirectories(outputDir);
                        break;
                    }
                case ARGUMENTS::TEST_MODE:
                    {
                        testMode    = true;
                        headerLines = {headerBaseString + " in test mode"};
                        std::cout << "Running in test mode\n";
                        break;
                    }
                case ARGUMENTS::OVERWRITE_EXISTING:
                    {
                        overwrite = true;
                    }
                default:
                    break;
            }
        }
        if (!(testMode || overwrite || outputDir == "."))
        {
            if (!codeUtils::directoryIsEmptyOrNonexistent(outputDir))
            {
                throw std::runtime_error("Output directory '" + outputDir +
                                         "' not empty. Please empty/remove it or run the program with the flag --" +
                                         overwriteFlag);
            }
        }
        // Convert command line inputs to legible variables
        std::ifstream infile(inputFile);
        std::string line;
        while (std::getline(infile, line))
        {
            try
            {
                if (line.empty())
                {
                    continue;
                }
                std::vector<std::string> splitLine = codeUtils::split(line, delimiter);
                if (splitLine.size() != 2)
                {
                    std::cerr << "Encountered problem when splitting this line >>>" << line << "<<< from file >>>"
                              << argv[1] << "<<< into an ID and carb string separated by your specified delimiter: >>>"
                              << argv[2] << "<<<\n";
                    std::exit(EXIT_FAILURE);
                }
                std::string inputSequence = splitLine.at(1);
                std::cout << "\n*********************\nBuilding " << inputSequence << "\n*********************\n";
                cdsCondensedSequence::carbohydrateBuilder carbBuilder(inputSequence);
                std::string inputGlycanID = splitLine.at(0);
                carbBuilder.GenerateSingle3DStructureDefaultFiles(outputDir, inputGlycanID);
            }
            catch (const std::runtime_error& error)
            {
                gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
                std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << error.what()
                          << std::endl;
            }
            catch (...)
            {
                gmml::log(__LINE__, __FILE__, gmml::ERR,
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
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
    catch (...)
    {
        std::string message =
            "carbohydrateBuilder: unexpected error. Please report how you got this to glycam@gmail.com";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
}
