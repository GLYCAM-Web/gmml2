#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/random.hpp"
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
        RNG_SEED,
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
        {ArgReq::optional,  ArgType::option,           RNG_SEED,    "rng-seed", ' ',            "value"},
        {ArgReq::optional,    ArgType::flag,               HELP,        "help", 'h',                 ""},
        {ArgReq::optional,    ArgType::flag,            VERSION,     "version", 'v',                 ""},
        {  ArgReq::hidden,    ArgType::flag,          TEST_MODE,   "test-mode", ' ',                 ""},
        {ArgReq::optional,    ArgType::flag, OVERWRITE_EXISTING, overwriteFlag, ' ',                 ""}
    };

    std::string programName = codeUtils::programName(argv);
    std::string baseDir     = codeUtils::toString(codeUtils::pathAboveCurrentExecutableDir());
    std::string SNFGDir     = codeUtils::SNFGSymbolsDir();

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
            std::cout << "Carbohydrate Builder & GMML2 version " << GMML_VERSION << "\n";
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
        std::string headerBaseString = "Produced by GMML2 (https://github.com/GLYCAM-Web/gmml2)";
        std::vector<std::string> headerLines {headerBaseString + " version " + std::string(GMML_VERSION)};
        uint64_t rngSeed = codeUtils::generateRandomSeed();
        bool testMode    = false;
        bool overwrite   = false;
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
                        break;
                    }
                case ARGUMENTS::OUTPUT_DIR:
                    {
                        outputDir = arg.value;
                        codeUtils::createDirectories(outputDir);
                        break;
                    }
                case ARGUMENTS::RNG_SEED:
                    {
                        std::optional<ulong> seed = codeUtils::parseUlong(arg.value);
                        if (!seed.has_value())
                        {
                            throw std::runtime_error("Error: rng-seed not a valid non-negative integer: " + arg.value);
                        }
                        rngSeed = seed.value();
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
                        break;
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

        std::string dotBaseDir = baseDir + "/";
        if (!testMode)
        {
            SNFGDir    = dotBaseDir + SNFGDir;
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
                std::vector<std::string> splitLine = codeUtils::split(line, delimiter);
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
        codeUtils::readFileLineByLine(inputFile, processLine);
        const cdsParameters::ParameterManager parameterManager = cdsParameters::loadParameters(baseDir);
        const codeUtils::SparseVector<double>& elementRadii    = MolecularMetadata::defaultVanDerWaalsRadii();
        pcg32 seedingRng(rngSeed);
        std::vector<uint64_t> rngSeeds = codeUtils::randomIntegers<uint64_t>(lines.size(), seedingRng);
        for (size_t n = 0; n < lines.size(); n++)
        {
            const SequenceInput& line = lines[n];
            pcg32 rng(rngSeeds[n]);
            try
            {
                std::cout << "\n*********************\nBuilding " << line.id << ": " << line.sequence
                          << "\n*********************\n";
                cdsCondensedSequence::SequenceData sequenceData = cdsCondensedSequence::parseSequence(line.sequence);
                graph::Graph traversable                        = graph::identity(sequenceData.graph);
                size_t nodeCount                                = traversable.nodes.indices.size();
                std::vector<size_t> parent(nodeCount, 0);
                std::vector<size_t> livingParent(nodeCount, 0);
                std::vector<bool> isBranchHead(nodeCount, false);
                std::vector<size_t> actualHead = codeUtils::indexVector(nodeCount);
                std::vector<size_t> actualTail = codeUtils::indexVector(nodeCount);
                std::vector<size_t> traversal  = {0};
                std::vector<bool> included(nodeCount, false);
                std::vector<bool> traversed(nodeCount, false);
                included[0]  = true;
                traversed[0] = true;
                for (size_t n = 0; n < traversal.size(); n++)
                {
                    size_t current                                 = traversal[n];
                    const cdsCondensedSequence::NodeType& nodeType = sequenceData.nodes.nodeType[current];
                    // means we're dealing with a probability node
                    if (std::holds_alternative<cdsCondensedSequence::ProbabilityNode>(nodeType))
                    {
                        const cdsCondensedSequence::ProbabilityNode& node =
                            std::get<cdsCondensedSequence::ProbabilityNode>(nodeType);
                        sequenceData.graph.nodeAlive[current] = false;
                        included[current] = codeUtils::uniformRandomDoubleWithinRange(rng, 0.0, 1.0) < node.probability;
                        size_t choice     = codeUtils::weightedRandomIndex(rng, node.weights);
                        size_t id         = node.heads[choice];
                        actualHead[current] = id;
                        actualTail[current] = node.tails[choice];
                        for (size_t k : node.heads)
                        {
                            parent[k]    = current;
                            traversed[k] = true;
                            included[k]  = false;
                        }
                        if (included[current])
                        {
                            traversal.push_back(id);
                            included[id] = true;
                        }
                    }
                    else if (std::holds_alternative<cdsCondensedSequence::BranchNode>(nodeType))
                    {
                        const cdsCondensedSequence::BranchNode& node =
                            std::get<cdsCondensedSequence::BranchNode>(nodeType);
                        sequenceData.graph.nodeAlive[current] = false;
                        size_t id                             = node.head;
                        parent[id]                            = current;
                        traversed[id]                         = true;
                        traversal.push_back(id);
                        included[id] = true;
                    }
                    const std::vector<size_t>& adjacencies = traversable.nodes.nodeAdjacencies[current];
                    for (size_t k = 0; k < adjacencies.size(); k++)
                    {
                        size_t neighbor = adjacencies[k];
                        if (!traversed[neighbor])
                        {
                            traversal.push_back(neighbor);
                            parent[neighbor]    = current;
                            traversed[neighbor] = true;
                            included[neighbor]  = true;
                        }
                    }
                }
                for (size_t n = traversal.size() - 1; n < traversal.size(); n--)
                {
                    size_t current      = traversal[n];
                    actualHead[current] = actualHead[actualHead[current]];
                    actualTail[current] = actualTail[actualTail[current]];
                }
                for (size_t n = 0; n < traversal.size(); n++)
                {
                    size_t current        = traversal[n];
                    livingParent[current] = included[parent[current]] ? parent[current] : livingParent[parent[current]];
                    const cdsCondensedSequence::NodeType& nodeType = sequenceData.nodes.nodeType[current];
                    if (std::holds_alternative<cdsCondensedSequence::ProbabilityNode>(nodeType))
                    {
                        // included[current] = included[current] && included[actualTail[parent[current]]];
                        // when dropped, residues that would connect to this in sequence will connect to parent instead
                        if (!included[current])
                        {
                            // livingParent[current] = livingParent[parent[current]];
                        }
                    }
                }
                sequenceData.graph.nodeAlive = codeUtils::vectorAnd(sequenceData.graph.nodeAlive, included);
                for (size_t n = 0; n < nodeCount; n++)
                {
                    std::cout << livingParent[n] << " " << parent[n] << " " << n << " (" << actualHead[n] << " - "
                              << actualTail[n] << ") " << included[n] << " " << sequenceData.nodes.fullString[n]
                              << "\n";
                }

                // add missing connections
                for (size_t nodeId : codeUtils::boolsToIndices(included))
                {
                    size_t n                                   = actualHead[nodeId];
                    size_t lastParent                          = livingParent[nodeId];
                    cdsCondensedSequence::NodeType& nodeType   = sequenceData.nodes.nodeType[nodeId];
                    cdsCondensedSequence::NodeType& parentType = sequenceData.nodes.nodeType[lastParent];
                    std::cout << parent[nodeId] << " " << nodeId << " ";
                    if (std::holds_alternative<cdsCondensedSequence::BranchNode>(parentType))
                    {
                        isBranchHead[n] = true;
                        cdsCondensedSequence::BranchNode& branch =
                            std::get<cdsCondensedSequence::BranchNode>(parentType);
                        std::string linkageName = sequenceData.nodes.linkage[n];
                        size_t dashPosition     = linkageName.find('-');
                        if (dashPosition == std::string::npos)
                        {
                            throw std::runtime_error("Missing dash in branch linkage name: " +
                                                     sequenceData.nodes.fullString[n]);
                        }
                        sequenceData.nodes.linkage[n] = linkageName.substr(0, dashPosition + 1) + branch.linkage;
                        graph::addEdge(sequenceData.graph, {actualTail[livingParent[lastParent]], n});
                        sequenceData.edges.names.push_back("");
                        std::cout << "branch " << actualTail[livingParent[lastParent]] << " " << n << " " << linkageName
                                  << " " << branch.linkage;
                    }
                    else if (std::holds_alternative<cdsCondensedSequence::ProbabilityNode>(nodeType) ||
                             !included[parent[nodeId]] ||
                             (actualTail[lastParent] != parent[nodeId]) && (actualHead[lastParent] != n))
                    {
                        graph::addEdge(sequenceData.graph, {actualTail[lastParent], n});
                        sequenceData.edges.names.push_back("");
                        std::cout << "skip " << actualTail[lastParent] << " " << n;
                    }
                    std::cout << "\n";
                }

                for (size_t n = 0; n < sequenceData.graph.edges.size(); n++)
                {
                    const std::array<size_t, 2>& edge = sequenceData.graph.edgeNodes[n];
                    size_t parent                     = edge[0];
                    size_t child                      = edge[1];
                    if (sequenceData.graph.nodeAlive[parent] && sequenceData.graph.nodeAlive[child])
                    {
                        bool isChildDeoxy = sequenceData.nodes.type[child] == cds::ResidueType::Deoxy;
                        bool isChildSugar = sequenceData.nodes.type[child] == cds::ResidueType::Sugar;
                        if (!isBranchHead[child] && !sequenceData.nodes.chainPosition[parent].empty())
                        {
                            std::string linkageName = sequenceData.nodes.linkage[child];
                            size_t dashPosition     = linkageName.find('-');
                            if (dashPosition != std::string::npos)
                            {
                                sequenceData.nodes.linkage[child] =
                                    linkageName.substr(0, dashPosition + 1) + sequenceData.nodes.chainPosition[parent];
                            }
                        }
                        // It remains internal if it's already been made internal, or if the child is not a deoxy.
                        sequenceData.nodes.isInternal[parent] = sequenceData.nodes.isInternal[parent] || !isChildDeoxy;
                        sequenceData.edges.names[n] =
                            sequenceData.nodes.configuration[child] + sequenceData.nodes.linkage[child];
                    }
                }
                for (size_t n = 0; n < nodeCount; n++)
                {
                    if (sequenceData.graph.nodeAlive[n])
                    {
                        std::cout << n << " " << sequenceData.nodes.fullString[n] << "\n";
                    }
                }
                for (size_t n = 0; n < sequenceData.graph.edgeNodes.size(); n++)
                {
                    auto& edge = sequenceData.graph.edgeNodes[n];
                    if (sequenceData.graph.nodeAlive[edge[0]] && sequenceData.graph.nodeAlive[edge[1]])
                    {
                        std::cout << edge[0] << " - " << edge[1] << " " << sequenceData.edges.names[n] << "\n";
                    }
                }
                cdsCondensedSequence::SequenceData pruned    = cdsCondensedSequence::pruned(sequenceData);
                cdsCondensedSequence::SequenceData reordered = cdsCondensedSequence::reordered(pruned);
                cdsCondensedSequence::Carbohydrate carbohydrate(parameterManager, elementRadii, reordered);
                carbohydrate.Generate3DStructureFiles(outputDir, line.id, headerLines);
                cdsCondensedSequence::GraphVizDotConfig config(dotBaseDir, SNFGDir, outputDir + "/" + line.id + ".dot");
                cdsCondensedSequence::printGraphViz(config, pruned);
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
