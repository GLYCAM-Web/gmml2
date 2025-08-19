#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/atomicConnectivity.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/glycoprotein/cdsInterface.hpp"
#include "include/glycoprotein/glycoproteinCreation.hpp"
#include "include/glycoprotein/sidechains.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/metadata/glycoprotein.hpp"
#include "include/metadata/sidechainRotamers.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/preprocess/pdbPreProcess.hpp"
#include "include/programs/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "include/programs/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/programs/GlycoproteinBuilder/overlapResolution.hpp"
#include "include/util/arguments.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
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
    using namespace gmml::gpbuilder;
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
        GlycoproteinBuilderInputs settings = readGPInputFile(inputFile);
        std::cout << "Reading input file complete, on to construction\n" << std::flush;

        const preprocess::ParameterManager parameterManager = preprocess::loadParameters(baseDir);
        pdb::PdbFile pdbFile = pdb::toPdbFile(settings.substrateFileName, pdb::modelsAsMolecules);
        if (settings.MDprep)
        {
            util::log(__LINE__, __FILE__, util::INF, "Performing MDPrep aka preprocessing.");
            preProcess(pdbFile, parameterManager, preprocess::defaultPreprocessorOptions);
        }

        util::SparseVector<double> elementRadii = vanDerWaalsRadii();
        const AminoAcidTable& aminoAcidTable = gmml::aminoAcidTable();
        const GlycosylationTable glycosylationTable = defaultGlycosylationTable();
        const DihedralAngleDataTable& dihedralAngleDataTable = gmml::dihedralAngleDataTable();
        const std::vector<Sphere> atomBounds =
            assembly::toAtomBounds(elementRadii, pdbFile.data.atoms.elements, pdbFile.data.atoms.coordinates);
        const assembly::Graph pdbGraph =
            assembly::createAssemblyGraph(pdbFile.data.assembly.indices, pdbFile.data.assembly.atomGraph);
        const assembly::Bounds bounds = assembly::toAssemblyBounds(pdbGraph, atomBounds);
        size_t glycoproteinAssemblyId = 0;
        Assembly* glycoprotein = getAssemblies(pdbFile)[glycoproteinAssemblyId];
        pdb::setIntraConnectivity(aminoAcidTable, pdbFile.data);
        pdb::setInterConnectivity(aminoAcidTable, pdbFile.data, bounds);
        util::log(__LINE__, __FILE__, util::INF, "Attaching Glycans To Glycosites.");
        std::vector<GlycosylationSite> glycosites =
            createGlycosites(pdbFile.data, glycoproteinAssemblyId, settings.glycositesInputVector);
        std::vector<Molecule*> glycans;
        std::vector<std::vector<ResidueLinkage>> glycosidicLinkages;
        addGlycansToProtein(
            glycans,
            glycosidicLinkages,
            glycosylationTable,
            parameterManager,
            elementRadii,
            dihedralAngleDataTable,
            pdbFile.data,
            glycoprotein,
            glycosites);
        util::log(__LINE__, __FILE__, util::INF, "Initialization of Glycoprotein builder complete!");

        std::vector<Molecule*> molecules = glycoprotein->getMolecules();
        GraphIndexData graphData = toIndexData(molecules);
        assembly::Graph graph = createCompleteAssemblyGraph(graphData);
        std::vector<size_t> glycanIndices = util::indicesOf(molecules, glycans);
        std::function<size_t(const GlycosylationSite&)> glycositeId = [](const GlycosylationSite& site)
        { return site.residueId; };
        std::vector<size_t> glycositeIds = util::vectorMap(glycositeId, glycosites);

        PartialLinkageData linkages = linkageData(dihedralAngleDataTable, graphData, graph, glycosidicLinkages);

        GlycoproteinAssembly assembly = addSidechainRotamers(
            aminoAcidTable,
            sidechainRotamers,
            toGlycoproteinAssemblyStructs(
                aminoAcidTable,
                elementRadii,
                pdbFile.data,
                graphData,
                graph,
                glycanIndices,
                glycositeIds,
                linkages,
                settings.ignoreHydrogen));
        if (settings.moveOverlappingSidechains)
        {
            assembly.data.atoms.includeInMainOverlapCheck = util::vectorAnd(
                assembly.data.atoms.includeInMainOverlapCheck,
                util::vectorNot(assembly.data.atoms.partOfMovableSidechain));
        }

        std::vector<bool> foundElements = gmml::foundElements(assembly.data.atoms.elements);

        OverlapSettings overlapSettings {
            potentialTable(elementRadii, foundElements),
            settings.overlapTolerance,
            settings.overlapRejectionThreshold,
        };

        ResolutionSettings resolution {
            settings.isDeterministic ? settings.seed : util::generateRandomSeed(),
            settings.persistCycles,
            settings.numberOfSamples,
            settings.useInitialGlycositeResidueConformation,
            settings.moveOverlappingSidechains,
            settings.deleteSitesUntilResolved,
            settings.MDprep,
            settings.MDprep};

        std::cout << "Resolving overlaps" << std::endl;
        std::vector<StructureStats> stats;
        resolveOverlaps(
            stats,
            sidechainRotamers,
            dihedralAngleDataTable,
            assembly,
            overlapSettings,
            resolution,
            outputDir,
            headerLines,
            numThreads);

        Summary summary = summarizeStats(assembly.graph, assembly.data, settings, resolution.rngSeed, stats);
        std::vector<util::TextVariant> textStructure = {
            util::TextVariant(util::TextHeader {2, "Glycoprotein Builder"}
             ),
            util::TextVariant(util::TextParagraph {headerLines}
             ),
            util::TextVariant(util::TextHeader {3, "Input"}
             ),
            util::TextVariant(
                util::TextParagraph {{"Filename: " + summary.filename, "Protein: " + summary.proteinFilename}}
             ),
            util::TextVariant(summary.parameterTable),
            util::TextVariant(util::TextHeader {3, "Structures"}
             ),
            util::TextVariant(summary.structuretable),
        };

        util::writeToFile(
            outputDir + "/summary.txt", [&](std::ostream& stream) { util::toTxt(stream, textStructure); });
        util::writeToFile(
            outputDir + "/summary.html", [&](std::ostream& stream) { util::toHtml(stream, textStructure); });
        util::writeToFile(
            outputDir + "/structures.csv",
            [&](std::ostream& stream) { util::toCsv(stream, ",", summary.structuretable); });
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
