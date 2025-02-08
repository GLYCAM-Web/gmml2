#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/threads.hpp"
#include "includes/version.h"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/sidechains.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <libgen.h>

int main(int argc, char* argv[])
{
    enum ARGUMENTS
    {
        INPUT_FILE,
        OUTPUT_DIR,
        HELP,
        VERSION,
        TEST_MODE,
        NUM_THREADS
    };

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,  INPUT_FILE,            "", ' ',       "input-file"},
        {ArgReq::optional, ArgType::unnamed,  OUTPUT_DIR,            "", ' ', "output-directory"},
        {ArgReq::optional,    ArgType::flag,        HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,     VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,   TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,  ArgType::option, NUM_THREADS, "num-threads", 'p',            "value"}
    };
    std::string programName = codeUtils::programName(argv);
    codeUtils::Path path    = codeUtils::toPath(argv[0]);
    size_t pathSize         = path.elements.size();
    if (pathSize < 3 || path.elements[pathSize - 2] != std::string("bin"))
    {
        throw std::runtime_error("expected application to be located in bin/" + path.elements[pathSize - 1] +
                                 ", actually located in " + programName);
    }

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
            std::cout << "Glycoprotein Builder & GMML version " << GMML_VERSION << "\n";
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
        std::string outputDir        = ".";
        std::string headerBaseString = "Produced by GMML (https://github.com/GLYCAM-Web/gmml2)";
        std::vector<std::string> headerLines {headerBaseString + " version " + std::string(GMML_VERSION)};
        int numThreads = 1;
        for (const auto& arg : arguments.args)
        {
            switch (arg.id)
            {
                case ARGUMENTS::INPUT_FILE:
                    {
                        inputFile = arg.value;
                        break;
                    }
                case ARGUMENTS::OUTPUT_DIR:
                    {
                        outputDir = arg.value;
                        struct stat info;
                        if (stat(outputDir.c_str(), &info) != 0)
                        {
                            std::cerr << "Folder " << outputDir << "/ does not exist and it isn't my job to make it.\n";
                            std::exit(EXIT_FAILURE);
                        }
                        break;
                    }
                case ARGUMENTS::TEST_MODE:
                    {
                        headerLines = {headerBaseString + " in test mode"};
                        std::cout << "Running in test mode\n";
                        break;
                    }
                case ARGUMENTS::NUM_THREADS:
                    {
                        std::optional<int> opt = codeUtils::parseInt(arg.value);
                        if (!opt.has_value() || opt.value() <= 0)
                        {
                            throw std::runtime_error(arg.value + " is not a valid value for " + arg.name +
                                                     ", must be a positive integer\n");
                        }
                        if (codeUtils::isOpenMpDefined())
                        {
                            numThreads = opt.value();
                            std::cout << "Number of threads set to " << numThreads << std::endl;
                        }
                        else
                        {
                            std::cout << "OpenMP not available, running on a single thread" << std::endl;
                        }
                        break;
                    }
                default:
                    break;
            }
        }
        MolecularMetadata::SidechainRotamerData sidechainRotamers;
        {
            std::string basePath = (path.absolute ? "/" : "") +
                                   codeUtils::join("/", codeUtils::take(path.elements.size() - 2, path.elements));
            std::string dunbrackLib = basePath + "/dat/dunbrack/sidechainRotamers.txt";
            sidechainRotamers       = MolecularMetadata::readSidechainRotamerData(dunbrackLib);
        }
        std::cout << "Input file is " << inputFile << "\n";
        glycoproteinBuilder::GlycoproteinBuilderInputs settings = glycoproteinBuilder::readGPInputFile(inputFile);
        std::cout << "Reading input file complete, on to construction\n" << std::flush;

        pdb::PdbFile pdbFile(settings.substrateFileName);
        if (!settings.skipMDPrep)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
            pdbFile.PreProcess(pdb::PreprocessorOptions());
        }

        cds::Assembly* glycoprotein                  = &pdbFile.mutableAssemblies().front();
        std::vector<cds::Residue*> gpInitialResidues = glycoprotein->getResidues();
        cds::setIntraConnectivity(gpInitialResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        std::vector<glycoproteinBuilder::GlycosylationSite> glycosites =
            glycoproteinBuilder::createGlycosites(glycoprotein, settings.glycositesInputVector);
        cds::setInterConnectivity(gpInitialResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");

        double selfWeight                                = 1000000.0;
        double proteinWeight                             = 1000.0;
        double glycanWeight                              = 1.0;
        glycoproteinBuilder::OverlapWeight overlapWeight = {proteinWeight, glycanWeight, selfWeight};

        std::vector<cds::Molecule*> molecules = glycoprotein->getMolecules();
        const glycoproteinBuilder::GlycoproteinAssembly initialAssembly =
            glycoproteinBuilder::toGlycoproteinAssemblyStructs(molecules, glycosites, overlapWeight);

        const glycoproteinBuilder::GlycoproteinAssembly assembly =
            settings.allowSidechainAdjustment ? addSidechainRotamers(sidechainRotamers, initialAssembly)
                                              : initialAssembly;

        std::cout << "Resolving overlaps" << std::endl;
        glycoproteinBuilder::resolveOverlaps(sidechainRotamers, overlapWeight, assembly, settings, outputDir,
                                             headerLines, numThreads);
    }
    catch (const std::runtime_error& error)
    {
        std::string errorMessage(error.what());
        std::string message = "glycoproteinBuilder: " + errorMessage;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        std::cout << message << "\n";
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
