#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <sstream>
#include <ostream>
#include <functional>

using cdsCondensedSequence::ParsedResidue;

namespace
{
    using cdsCondensedSequence::ResidueData;
    using cdsCondensedSequence::SequenceData;
    using cdsCondensedSequence::SequencePrintConfig;
    using codeUtils::squareBrackets;

    auto residueNames = [](const SequenceData& sequence)
    {
        const ResidueData& residues = sequence.residues;
        size_t residueCount         = residues.name.size();
        std::vector<std::string> result;
        result.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            result.push_back(residues.isomer[n] + residues.name[n] + residues.ringType[n] + residues.modifier[n] +
                             residues.ringShape[n]);
        }
        return result;
    };

    auto iupacNames = [](const SequenceData& sequence)
    {
        const ResidueData& residues = sequence.residues;
        size_t residueCount         = residues.name.size();
        std::vector<std::string> result;
        result.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            cds::ResidueType type = residues.type[n];
            // Iupac doesn't write -OH, but does include reducing terminal when it's sugar-sugar ano-ano
            if (type == cds::ResidueType::Sugar)
            {
                result.push_back(residues.name[n] + residues.modifier[n]);
            }
            else if (codeUtils::contains({cds::ResidueType::Deoxy, cds::ResidueType::Derivative}, type))
            {
                result.push_back(residues.isomer[n] + residues.name[n] + residues.ringType[n] + residues.modifier[n] +
                                 residues.ringShape[n]);
            }
            else
            {
                result.push_back("");
            }
        }
        return result;
    };

    auto noResidueLabels = [](const SequenceData& sequence)
    {
        return std::vector<std::string>(sequence.graph.nodes.size(), "");
    };

    auto residueLabels = [](const SequenceData& sequence)
    {
        const ResidueData& residues = sequence.residues;
        size_t residueCount         = residues.name.size();
        std::vector<std::string> result;
        result.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            result.push_back("&Label=residue-" + std::to_string(n + 1) + ";");
        }
        return result;
    };

    auto residueConfigurations = [](const SequenceData& sequence)
    {
        return sequence.residues.configuration;
    };

    auto noResidueConfigurations = [](const SequenceData& sequence)
    {
        return std::vector<std::string>(sequence.graph.nodes.size(), "");
    };

    auto edgeNames = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>&)
    {
        return sequence.edges.names;
    };

    auto noEdgeLabels = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>&)
    {
        return std::vector<std::string>(sequence.graph.edges.size(), "");
    };

    auto edgeLabels = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>& nodeEdges)
    {
        size_t edgeCount = sequence.graph.edges.size();
        size_t index     = 0;
        std::vector<size_t> edgeIndices(edgeCount, 0);
        for (size_t nodeId : sequence.graph.nodes)
        {
            for (size_t edgeId : nodeEdges[nodeId])
            {
                edgeIndices[edgeId] = index;
                index++;
            }
        }

        std::vector<std::string> result;
        result.reserve(edgeCount);
        for (size_t n : sequence.graph.edges)
        {
            result.push_back("&Label=link-" + std::to_string(edgeIndices[n]) + ";");
        }
        return result;
    };

    auto noSort = [](const SequenceData&, const std::vector<size_t>& edgeIds)
    {
        return edgeIds;
    };

    auto sortEdges = [](const SequenceData& sequence, const std::vector<size_t>& edgeIds)
    {
        return cdsCondensedSequence::edgesSortedByLink(sequence, edgeIds);
    };

    auto noLinkageBrackets = [](const SequenceData& sequence)
    {
        return std::vector<codeUtils::Brackets>(sequence.graph.nodes.size(), codeUtils::noBrackets);
    };

    auto iupacLinkageBrackets = [](const SequenceData& sequence)
    {
        std::vector<codeUtils::Brackets> result(sequence.graph.nodes.size(), {"", ""});
        for (size_t nodeId : sequence.graph.nodes)
        {
            if (sequence.residues.type[nodeId] != cds::ResidueType::Aglycone)
            {
                result[nodeId].open = "(";
            }
        }
        for (size_t edgeId : sequence.graph.edges)
        {
            const std::array<size_t, 2>& nodes = sequence.graph.edgeNodes[edgeId];
            if (!codeUtils::contains({sequence.residues.type[nodes[0]], sequence.residues.type[nodes[1]]},
                                     cds::ResidueType::Aglycone))
            {
                result[nodes[1]].close = ")";
            }
        }
        return result;
    };

    struct SequencePrintData
    {
        codeUtils::Brackets derivativeBrackets;
        std::string derivativeSeparator;
        std::vector<std::vector<size_t>> residueChildEdges;
        std::vector<cds::ResidueType> residueTypes;
        std::vector<std::string> residueNames;
        std::vector<std::string> parentlessResidueLabel;
        std::vector<codeUtils::Brackets> residueLinkageBrackets;
        std::vector<std::string> edgeNames;
    };

    std::vector<size_t> childEdges(const SequenceData& sequence, size_t residueId)
    {
        size_t edgeCount = sequence.graph.edges.size();
        std::vector<size_t> result;
        result.reserve(edgeCount);
        for (size_t n = 0; n < edgeCount; n++)
        {
            if (sequence.graph.edgeNodes[n][0] == residueId)
            {
                result.push_back(n);
            }
        }
        return result;
    }

    std::vector<size_t> parentEdges(const SequenceData& sequence, size_t residueId)
    {
        size_t edgeCount = sequence.graph.edges.size();
        std::vector<size_t> result;
        result.reserve(edgeCount);
        for (size_t n = 0; n < edgeCount; n++)
        {
            if (sequence.graph.edgeNodes[n][1] == residueId)
            {
                result.push_back(n);
            }
        }
        return result;
    }

    int recurvePrint(const SequencePrintData& data, const graph::Database& graph, size_t residueId, size_t parentEdgeId,
                     int branchStackSize, std::vector<std::string>& output)
    {
        const std::vector<size_t>& edges = data.residueChildEdges[residueId];
        size_t numberOfNeighbors         = edges.size();
        // Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
        std::string outputResidueString  = data.residueNames[residueId];
        std::vector<std::string> derivatives;
        for (size_t edgeId : edges)
        {
            size_t neighbor = graph.edgeNodes[edgeId][1];
            if (codeUtils::contains({cds::ResidueType::Derivative, cds::ResidueType::Deoxy},
                                    data.residueTypes[neighbor]))
            {
                --numberOfNeighbors;
                derivatives.push_back(data.edgeNames[edgeId] + data.residueNames[neighbor]);
            }
        }
        if (!derivatives.empty())
        {
            outputResidueString += data.derivativeBrackets.open;
            outputResidueString += codeUtils::join(
                data.derivativeSeparator, codeUtils::reverse(derivatives)); // order should be 2S,6S, not 6S,2S.
            outputResidueString += data.derivativeBrackets.close;
        }
        outputResidueString += data.residueLinkageBrackets[residueId].open;
        bool hasParentEdge  = parentEdgeId < graph.edges.size();
        outputResidueString += (hasParentEdge ? data.edgeNames[parentEdgeId] : data.parentlessResidueLabel[residueId]);
        outputResidueString += data.residueLinkageBrackets[residueId].close;
        output.push_back(outputResidueString);
        // End of a branch check
        if (numberOfNeighbors == 0 && branchStackSize > 0)
        {
            output.push_back(squareBrackets.open);
            --branchStackSize;
        }
        size_t loopCount = 0;
        for (size_t edgeId : edges)
        {
            size_t neighbor = graph.edgeNodes[edgeId][1];
            if (!codeUtils::contains({cds::ResidueType::Derivative, cds::ResidueType::Deoxy},
                                     data.residueTypes[neighbor]))
            {
                ++loopCount;
                if (loopCount < numberOfNeighbors)
                {
                    output.push_back(squareBrackets.close);
                    ++branchStackSize;
                }
                branchStackSize = recurvePrint(data, graph, neighbor, edgeId, branchStackSize, output);
            }
        }
        return branchStackSize;
    }
} // namespace

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    auto residueLink = [&](size_t n)
    {
        return cdsCondensedSequence::getLink(sequence.residues.type[n], sequence.residues.linkage[n]);
    };
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        return residueLink(sequence.graph.edgeNodes[n][1]) > residueLink(sequence.graph.edgeNodes[k][1]);
    };

    return codeUtils::sortedBy(compare, edgeIds);
}

std::string cdsCondensedSequence::printGraphViz(GraphVizDotConfig& configs, const SequenceData& sequence)
{
    std::vector<std::string> derivatives         = sequenceDerivatives(sequence);
    std::vector<std::string> monosaccharideNames = sequenceMonosaccharideNames(sequence);
    graph::Graph graph                           = condensedSequenceGraph(sequence);

    std::vector<GraphVizResidueNode> nodes;
    nodes.reserve(nodeCount(graph));
    std::vector<GraphVizLinkage> linkages;
    linkages.reserve(edgeCount(graph));

    for (size_t n = 0; n < nodeCount(graph); n++)
    {
        size_t index                    = sourceNodeIndex(graph, n);
        std::string& monosaccharideName = monosaccharideNames[index];
        GraphVizImage image =
            findImage(configs, monosaccharideName, (sequence.residues.ringType[index] == "f") ? "f" : "");
        nodes.push_back({n, image, monosaccharideName, derivatives[index]});
    }

    for (size_t n = 0; n < edgeCount(graph); n++)
    {
        std::array<size_t, 2>& adj = graph.edges.nodeAdjacencies[n];
        size_t parent              = adj[0];
        size_t child               = adj[1];
        size_t childIndex          = sourceNodeIndex(graph, child);
        std::string label          = (configs.show_config_labels_ ? sequence.residues.configuration[childIndex] : "") +
                            (configs.show_position_labels_ ? sequence.residues.linkage[childIndex] : "");
        linkages.push_back({
            {child, parent},
            label
        });
    }

    std::stringstream ss;
    ss << "graph G {graph [splines=false dpi=" << configs.dpi_ << " outputorder=\"edgesfirst\"];\n";
    ss << "node [shape=\"none\" fontname=DejaVuSans labelfontsize=12 ";
    ss << "label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];\n";
    ss << "edge [labelfontsize=12 fontname=DejaVuSans labeldistance=1.2 labelangle=320.0];\n";
    ss << "rankdir=LR nodesep=\"0.05\" ranksep=\"0.8\";\n";
    for (size_t n = 0; n < nodes.size(); n++)
    {
        size_t nodeIndex = sourceNodeIndex(graph, n);
        bool isAglycone  = sequence.residues.type[nodeIndex] == cds::ResidueType::Aglycone;
        ss << (isAglycone ? graphVizAglyconeNode(nodes[n]) : graphVizSugarNode(nodes[n]));
    }
    for (auto& linkage : linkages)
    {
        ss << graphVizLinkageLine(nodes, linkage);
    }
    ss << "}\n";
    std::string str = ss.str();
    // Open and overwrite.
    codeUtils::writeToFile(configs.file_name_,
                           [&str](std::ostream& stream)
                           {
                               stream << str;
                           });
    return str;
}

std::vector<ParsedResidue*>
cdsCondensedSequence::parsedResiduesOrderedByConnectivity(std::vector<ParsedResidue*> residues)
{
    std::vector<ParsedResidue*> rawResidues;
    // Go via Graph so order decided by connectivity, depth first traversal:
    glygraph::Graph<cds::Residue> sequenceGraph(terminalResidue(residues));
    for (auto& node : sequenceGraph.getNodes())
    {
        rawResidues.push_back(codeUtils::erratic_cast<ParsedResidue*>(node->getDerivedClass()));
    }
    return rawResidues;
}

void cdsCondensedSequence::setIndexByConnectivity(std::vector<ParsedResidue*> residues)
{
    unsigned long long linkIndex    = 0; // Convention to start form 0 for linkages.
    unsigned long long residueIndex = 1; // Convention to start from 1 for residues.
    for (auto& residue : parsedResiduesOrderedByConnectivity(residues))
    {
        residue->setIndex(residueIndex);
        residue->setNumber(residueIndex); // ToDo temporary, switch to using number here. Keep index as a gmml internal
                                          // thing, never shown to user.
        ++residueIndex;
        for (auto& edge : residue->getInEdges())
        {
            edge->setIndex(linkIndex);
            ++linkIndex;
        }
    }
    return;
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::defaultConfig()
{
    return {sortEdges,
            codeUtils::squareBrackets,
            ",",
            residueNames,
            noResidueLabels,
            residueConfigurations,
            noLinkageBrackets,
            edgeNames,
            noEdgeLabels};
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::defaultConfigUnsorted()
{
    SequencePrintConfig config = defaultConfig();
    config.edgeOrder           = noSort;
    return config;
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::defaultConfigLabelled()
{
    SequencePrintConfig config     = defaultConfig();
    config.parentlessResidueLabels = noResidueConfigurations;
    config.residueLabels           = residueLabels;
    config.edgeLabels              = edgeLabels;
    return config;
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::iupacConfig()
{
    return {
        sortEdges, codeUtils::noBrackets, "", iupacNames, noResidueLabels, residueConfigurations, iupacLinkageBrackets,
        edgeNames, noEdgeLabels};
}

std::string cdsCondensedSequence::printSequence(const SequenceData& sequence, const SequencePrintConfig& config)
{
    auto joinStrings = [](const std::vector<std::string>& avec, const std::vector<std::string>& bvec)
    {
        std::vector<std::string> result;
        result.reserve(avec.size());
        for (size_t n = 0; n < avec.size(); n++)
        {
            result.push_back(avec[n] + bvec[n]);
        }
        return result;
    };
    size_t noEdgeId = size_t(-1);
    std::vector<std::vector<size_t>> residueChildEdges;
    std::vector<std::vector<size_t>> residueParentEdges;
    for (size_t nodeId : sequence.graph.nodes)
    {
        residueChildEdges.push_back(config.edgeOrder(sequence, childEdges(sequence, nodeId)));
        residueParentEdges.push_back(parentEdges(sequence, nodeId));
    }
    std::vector<std::string> residueNames = joinStrings(config.residueNames(sequence), config.residueLabels(sequence));
    std::vector<std::string> edgeNames =
        joinStrings(config.edgeNames(sequence, residueParentEdges), config.edgeLabels(sequence, residueParentEdges));
    SequencePrintData data = {config.derivativeBrackets,
                              config.derivativeSeparator,
                              residueChildEdges,
                              sequence.residues.type,
                              residueNames,
                              config.parentlessResidueLabels(sequence),
                              config.residueLinkageBrackets(sequence),
                              edgeNames};
    std::vector<std::string> output;
    recurvePrint(data, sequence.graph, 0, noEdgeId, 0, output);
    std::reverse(output.begin(), output.end()); // Reverse order, as it starts from terminal.
    std::stringstream ss;
    for (auto& label : output)
    {
        ss << label;
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::reordered(const SequenceData& sequence)
{
    size_t residueCount           = sequence.residues.name.size();
    std::vector<size_t> edgeOrder = edgesSortedByLink(sequence, codeUtils::indexVector(sequence.graph.edges));
    std::vector<std::array<size_t, 2>> reorderedEdges = codeUtils::indicesToValues(sequence.graph.edgeNodes, edgeOrder);

    size_t current                   = 0;
    std::vector<size_t> residueOrder = {current};
    residueOrder.reserve(residueCount);
    std::vector<bool> traversed = codeUtils::indicesToBools(residueCount, residueOrder);
    auto fromNode               = [&](size_t node)
    {
        std::vector<size_t> result;
        result.reserve(reorderedEdges.size());
        for (auto& edge : sequence.graph.edgeNodes)
        {
            if (edge[1] == node)
            {
                result.push_back(edge[0]);
            }
        }
        for (auto& edge : reorderedEdges)
        {
            if (edge[0] == node)
            {
                result.push_back(edge[1]);
            }
        }
        return codeUtils::reverse(result);
    };
    std::vector<size_t> nodesToTraverse;
    codeUtils::insertInto(nodesToTraverse, fromNode(current));
    while (!nodesToTraverse.empty())
    {
        current = nodesToTraverse.back();
        nodesToTraverse.pop_back();
        if (!traversed[current])
        {
            residueOrder.push_back(current);
            codeUtils::insertInto(nodesToTraverse, fromNode(current));
            traversed[current] = true;
        }
    }

    std::vector<size_t> missed = codeUtils::boolsToIndices(codeUtils::vectorNot(traversed));
    if (missed.size() > 0)
    {
        throw std::runtime_error("Error: sequence graph not fully connected");
    }

    std::vector<size_t> invertedResidueOrder(residueCount, -1);
    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        invertedResidueOrder[residueOrder[n]] = n;
    }

    graph::Database resultGraph;

    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        graph::addNode(resultGraph);
    }

    std::vector<size_t> edgeIndexOrder = codeUtils::indexVector(sequence.graph.edges);

    for (size_t n : edgeIndexOrder)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    ResidueData residues = {codeUtils::indicesToValues(sequence.residues.fullString, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.type, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.name, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.linkage, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringType, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.configuration, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isomer, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringShape, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.modifier, residueOrder)};

    return SequenceData {resultGraph, residues, sequence.edges};
}

void cdsCondensedSequence::createParsedResidues(
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs, const SequenceData& sequence)
{
    size_t residueCount = sequence.residues.name.size();
    residuePtrs.reserve(residueCount);
    for (size_t n = 0; n < residueCount; n++)
    {
        residuePtrs.emplace_back(
            std::make_unique<cdsCondensedSequence::ParsedResidue>(cdsCondensedSequence::ParsedResidueComponents {
                sequence.residues.fullString[n],
                sequence.residues.type[n],
                sequence.residues.name[n],
                sequence.residues.linkage[n],
                sequence.residues.ringType[n],
                sequence.residues.configuration[n],
                sequence.residues.isomer[n],
                {sequence.residues.ringShape[n], sequence.residues.modifier[n]}
        }));
    }
    for (size_t n = 0; n < sequence.graph.edgeNodes.size(); n++)
    {
        auto& edge = sequence.graph.edgeNodes[n];
        residuePtrs[edge[1]].get()->addParent(sequence.edges.names[n], residuePtrs[edge[0]].get());
    }
}

void cdsCondensedSequence::sortResidueEdges(
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs)
{
    for (auto& residue : residuePtrs)
    {
        residue.get()->sortOutEdgesBySourceTObjectComparator();
    }
}

ParsedResidue* cdsCondensedSequence::terminalResidue(std::vector<ParsedResidue*> residues)
{
    return residues.front();
}
