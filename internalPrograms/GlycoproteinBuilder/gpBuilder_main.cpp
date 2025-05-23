#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/cdsInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/sidechains.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/overlapResolution.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreProcess.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/threads.hpp"
#include "includes/version.h"

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
        NUM_THREADS,
        OVERWRITE_EXISTING
    };

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    const std::string overwriteFlag                    = "overwrite-existing-files";
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,         INPUT_FILE,            "", ' ',       "input-file"},
        {ArgReq::optional, ArgType::unnamed,         OUTPUT_DIR,            "", ' ', "output-directory"},
        {ArgReq::optional,    ArgType::flag,               HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,            VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,          TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,  ArgType::option,        NUM_THREADS, "num-threads", 't',            "value"},
        {ArgReq::optional,    ArgType::flag, OVERWRITE_EXISTING, overwriteFlag, ' ',                 ""}
    };
    std::string programName = codeUtils::programName(argv);
    std::string baseDir     = codeUtils::toString(codeUtils::pathAboveCurrentExecutableDir());

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
            std::cout << "Glycoprotein Builder & GMML2 version " << GMML_VERSION << "\n";
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
        std::string headerBaseString = "Produced by GMML2 (https://github.com/GLYCAM-Web/gmml2)";
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

        MolecularMetadata::SidechainRotamerData sidechainRotamers;
        {
            std::string dunbrackLib = baseDir + "/dat/dunbrack/sidechainRotamers.txt";
            sidechainRotamers       = MolecularMetadata::readSidechainRotamerData(dunbrackLib);
        }
        std::cout << "Input file is " << inputFile << "\n";
        glycoproteinBuilder::GlycoproteinBuilderInputs settings = glycoproteinBuilder::readGPInputFile(inputFile);
        std::cout << "Reading input file complete, on to construction\n" << std::flush;

        const cdsParameters::ParameterManager parameterManager = cdsParameters::loadParameters(baseDir);
        pdb::PdbFile pdbFile(settings.substrateFileName, {pdb::InputType::modelsAsMolecules, false});
        if (settings.MDprep)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
            pdbFile.PreProcess(parameterManager, pdb::PreprocessorOptions());
        }

        codeUtils::SparseVector<double> elementRadii                         = MolecularMetadata::vanDerWaalsRadii();
        const MolecularMetadata::AminoAcidTable& aminoAcidTable              = MolecularMetadata::aminoAcidTable();
        const GlycamMetadata::DihedralAngleDataTable& dihedralAngleDataTable = GlycamMetadata::dihedralAngleDataTable();
        cds::Assembly* glycoprotein                                          = pdbFile.getAssemblies().front();
        pdb::setIntraConnectivity(aminoAcidTable, pdbFile.data);
        pdb::setInterConnectivity(aminoAcidTable, pdbFile.data);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        std::vector<glycoproteinBuilder::GlycosylationSite> glycosites =
            glycoproteinBuilder::createGlycosites(pdbFile.data, glycoprotein, settings.glycositesInputVector);
        std::vector<cdsCondensedSequence::Carbohydrate*> glycans = glycoproteinBuilder::addGlycansToProtein(
            parameterManager, elementRadii, dihedralAngleDataTable, pdbFile.data, glycoprotein, glycosites);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");

        std::vector<cds::Molecule*> molecules = glycoprotein->getMolecules();

        glycoproteinBuilder::GlycoproteinAssembly assembly = addSidechainRotamers(
            aminoAcidTable, sidechainRotamers,
            glycoproteinBuilder::toGlycoproteinAssemblyStructs(
                aminoAcidTable, dihedralAngleDataTable, elementRadii, pdbFile.data, molecules, glycosites, glycans,
                settings.overlapTolerance, settings.overlapRejectionThreshold, settings.ignoreHydrogen));
        if (settings.moveOverlappingSidechains)
        {
            assembly.data.atoms.includeInMainOverlapCheck =
                codeUtils::vectorAnd(assembly.data.atoms.includeInMainOverlapCheck,
                                     codeUtils::vectorNot(assembly.data.atoms.partOfMovableSidechain));
        }

        std::cout << "Resolving overlaps" << std::endl;
        uint64_t rngSeed = settings.isDeterministic ? settings.seed : codeUtils::generateRandomSeed();
        glycoproteinBuilder::resolveOverlaps(sidechainRotamers, assembly, settings, rngSeed, outputDir, headerLines,
                                             numThreads);
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
