#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/metadata/sidechainRotamers.hpp"
#include "include/pdb/atomicConnectivity.hpp"
#include "include/pdb/pdbFile.hpp"
#include "include/pdb/pdbPreProcess.hpp"
#include "include/programs/GlycoproteinBuilder/cdsInterface.hpp"
#include "include/programs/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "include/programs/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/programs/GlycoproteinBuilder/overlapResolution.hpp"
#include "include/programs/GlycoproteinBuilder/sidechains.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/parsing.hpp"
#include "include/util/random.hpp"
#include "include/util/strings.hpp"
#include "include/util/threads.hpp"
#include "include/version.h"

#include <iostream>
#include <libgen.h>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    enum ARGUMENTS
    {
        INPUT_FILE,
        OUTPUT_DIR,
        HELP,
        VERSION,
        TEST_MODE,
        NUM_THREADS,
        OVERWRITE_EXISTING
    };

    using namespace gmml;
    using util::ArgReq;
    using util::ArgType;
    const std::string overwriteFlag = "overwrite-existing-files";
    std::vector<util::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,         INPUT_FILE,            "", ' ',       "input-file"},
        {ArgReq::optional, ArgType::unnamed,         OUTPUT_DIR,            "", ' ', "output-directory"},
        {ArgReq::optional,    ArgType::flag,               HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,            VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,          TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,  ArgType::option,        NUM_THREADS, "num-threads", 't',            "value"},
        {ArgReq::optional,    ArgType::flag, OVERWRITE_EXISTING, overwriteFlag, ' ',                 ""}
    };
    std::string programName = util::programName(argv);
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());

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
            std::cout << "Glycoprotein Builder & GMML2 version " << GMML_VERSION << "\n";
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
        std::string outputDir = ".";
        std::string headerBaseString = "Produced by GMML2 (https://github.com/GLYCAM-Web/gmml2)";
        std::vector<std::string> headerLines {headerBaseString + " version " + std::string(GMML_VERSION)};
        int numThreads = 1;
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
                case ARGUMENTS::NUM_THREADS:
                    {
                        std::optional<int> opt = util::parseInt(arg.value);
                        if (!opt.has_value() || opt.value() <= 0)
                        {
                            throw std::runtime_error(
                                arg.value + " is not a valid value for " + arg.name + ", must be a positive integer\n");
                        }
                        if (util::isOpenMpDefined())
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
            if (!util::directoryIsEmptyOrNonexistent(outputDir))
            {
                throw std::runtime_error(
                    "Output directory '" + outputDir +
                    "' not empty. Please empty/remove it or run the program with the flag --" + overwriteFlag);
            }
        }

        SidechainRotamerData sidechainRotamers;
        {
            std::string dunbrackLib = baseDir + "/dat/dunbrack/sidechainRotamers.txt";
            sidechainRotamers = readSidechainRotamerData(dunbrackLib);
        }
        std::cout << "Input file is " << inputFile << "\n";
        gpbuilder::GlycoproteinBuilderInputs settings = gpbuilder::readGPInputFile(inputFile);
        std::cout << "Reading input file complete, on to construction\n" << std::flush;

        const ParameterManager parameterManager = loadParameters(baseDir);
        pdb::PdbFile pdbFile = pdb::toPdbFile(settings.substrateFileName, pdb::modelsAsMolecules);
        if (settings.MDprep)
        {
            util::log(__LINE__, __FILE__, util::INF, "Performing MDPrep aka preprocessing.");
            preProcess(pdbFile, parameterManager, pdb::defaultPreprocessorOptions);
        }

        util::SparseVector<double> elementRadii = vanDerWaalsRadii();
        const AminoAcidTable& aminoAcidTable = gmml::aminoAcidTable();
        const DihedralAngleDataTable& dihedralAngleDataTable = gmml::dihedralAngleDataTable();
        Assembly* glycoprotein = getAssemblies(pdbFile).front();
        pdb::setIntraConnectivity(aminoAcidTable, pdbFile.data);
        pdb::setInterConnectivity(aminoAcidTable, pdbFile.data);
        util::log(__LINE__, __FILE__, util::INF, "Attaching Glycans To Glycosites.");
        std::vector<gpbuilder::GlycosylationSite> glycosites =
            createGlycosites(pdbFile.data, glycoprotein, settings.glycositesInputVector);
        std::vector<Carbohydrate*> glycans = addGlycansToProtein(
            parameterManager, elementRadii, dihedralAngleDataTable, pdbFile.data, glycoprotein, glycosites);
        util::log(__LINE__, __FILE__, util::INF, "Initialization of Glycoprotein builder complete!");

        std::vector<Molecule*> molecules = glycoprotein->getMolecules();

        gpbuilder::GlycoproteinAssembly assembly = addSidechainRotamers(
            aminoAcidTable,
            sidechainRotamers,
            toGlycoproteinAssemblyStructs(
                aminoAcidTable,
                dihedralAngleDataTable,
                elementRadii,
                pdbFile.data,
                molecules,
                glycosites,
                glycans,
                settings.overlapTolerance,
                settings.overlapRejectionThreshold,
                settings.ignoreHydrogen));
        if (settings.moveOverlappingSidechains)
        {
            assembly.data.atoms.includeInMainOverlapCheck = util::vectorAnd(
                assembly.data.atoms.includeInMainOverlapCheck,
                util::vectorNot(assembly.data.atoms.partOfMovableSidechain));
        }

        std::cout << "Resolving overlaps" << std::endl;
        uint64_t rngSeed = settings.isDeterministic ? settings.seed : util::generateRandomSeed();
        resolveOverlaps(sidechainRotamers, assembly, settings, rngSeed, outputDir, headerLines, numThreads);
    }
    catch (const std::runtime_error& error)
    {
        std::string errorMessage(error.what());
        std::string message = "glycoproteinBuilder: " + errorMessage;
        util::log(__LINE__, __FILE__, util::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
    catch (...)
    {
        std::string message =
            "glycoproteinBuilder: unexpected error. Please report how you got this to glycam@gmail.com";
        util::log(__LINE__, __FILE__, util::ERR, message);
        std::cout << message << "\n";
        std::exit(1);
    }
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
