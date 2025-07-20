#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"

#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceUtil.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"

#include <functional>
#include <ostream>
#include <sstream>

namespace
{
    using cdsCondensedSequence::ResidueData;
    using cdsCondensedSequence::SequenceData;
    using cdsCondensedSequence::SequencePrintConfig;
    using codeUtils::squareBrackets;

    auto residueNames = [](const SequenceData& sequence)
    {
        size_t residueCount = nodeCount(sequence.graph);
        std::vector<std::string> result;
        result.reserve(residueCount);
        const ResidueData& residues = sequence.residues;
        for (size_t n = 0; n < residueCount; n++)
        {
            result.push_back(
                residues.preIsomerModifier[n] + residues.isomer[n] + residues.name[n] + residues.ringType[n] +
                residues.modifier[n] + residues.ringShape[n]);
        }
        return result;
    };

    auto iupacNames = [](const SequenceData& sequence)
    {
        size_t residueCount = nodeCount(sequence.graph);
        std::vector<std::string> result;
        result.reserve(residueCount);
        const ResidueData& residues = sequence.residues;
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
                result.push_back(
                    residues.preIsomerModifier[n] + residues.isomer[n] + residues.name[n] + residues.ringType[n] +
                    residues.modifier[n] + residues.ringShape[n]);
            }
            else
            {
                result.push_back("");
            }
        }
        return result;
    };

    auto noResidueLabels = [](const SequenceData& sequence)
    { return std::vector<std::string>(nodeCount(sequence.graph), ""); };

    auto residueLabels = [](const SequenceData& sequence)
    {
        size_t residueCount = nodeCount(sequence.graph);
        std::vector<std::string> result;
        result.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            result.push_back("&Label=residue-" + std::to_string(n + 1) + ";");
        }
        return result;
    };

    auto residueConfigurations = [](const SequenceData& sequence) { return sequence.residues.configuration; };

    auto noResidueConfigurations = [](const SequenceData& sequence)
    { return std::vector<std::string>(nodeCount(sequence.graph), ""); };

    auto edgeNames = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>&)
    { return sequence.edges.names; };

    auto noEdgeLabels = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>&)
    { return std::vector<std::string>(edgeCount(sequence.graph), ""); };

    auto edgeLabels = [](const SequenceData& sequence, const std::vector<std::vector<size_t>>& nodeEdges)
    {
        size_t edgeCount = graph::edgeCount(sequence.graph);
        size_t index = 0;
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

    auto noSort = [](const SequenceData&, const std::vector<size_t>& edgeIds) { return edgeIds; };

    auto sortEdges = [](const SequenceData& sequence, const std::vector<size_t>& edgeIds)
    { return cdsCondensedSequence::edgesSortedByLink(sequence, edgeIds); };

    auto noLinkageBrackets = [](const SequenceData& sequence)
    { return std::vector<codeUtils::Brackets>(nodeCount(sequence.graph), codeUtils::noBrackets); };

    auto iupacLinkageBrackets = [](const SequenceData& sequence)
    {
        std::vector<codeUtils::Brackets> result(nodeCount(sequence.graph), {"", ""});
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
            if (!codeUtils::contains(
                    {sequence.residues.type[nodes[0]], sequence.residues.type[nodes[1]]}, cds::ResidueType::Aglycone))
            {
                result[nodes[1]].close = ")";
            }
        }
        return result;
    };

    std::vector<size_t> childEdges(const SequenceData& sequence, size_t residueId)
    {
        size_t edgeCount = graph::edgeCount(sequence.graph);
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
        size_t edgeCount = graph::edgeCount(sequence.graph);
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

    int recurvePrint(
        const SequencePrintData& data,
        const graph::Database& graph,
        size_t residueId,
        size_t parentEdgeId,
        int branchStackSize,
        std::vector<std::string>& output)
    {
        const std::vector<size_t>& edges = data.residueChildEdges[residueId];
        size_t numberOfNeighbors = edges.size();
        // Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
        std::string outputResidueString = data.residueNames[residueId];
        std::vector<std::string> derivatives;
        for (size_t edgeId : edges)
        {
            size_t neighbor = graph.edgeNodes[edgeId][1];
            if (codeUtils::contains(
                    {cds::ResidueType::Derivative, cds::ResidueType::Deoxy}, data.residueTypes[neighbor]))
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
        bool hasParentEdge = parentEdgeId < edgeCount(graph);
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
            if (!codeUtils::contains(
                    {cds::ResidueType::Derivative, cds::ResidueType::Deoxy}, data.residueTypes[neighbor]))
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

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::defaultConfig()
{
    return {
        sortEdges,
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
    config.edgeOrder = noSort;
    return config;
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::defaultConfigLabelled()
{
    SequencePrintConfig config = defaultConfig();
    config.parentlessResidueLabels = noResidueConfigurations;
    config.residueLabels = residueLabels;
    config.edgeLabels = edgeLabels;
    return config;
}

cdsCondensedSequence::SequencePrintConfig cdsCondensedSequence::iupacConfig()
{
    return {
        sortEdges,
        codeUtils::noBrackets,
        "",
        iupacNames,
        noResidueLabels,
        residueConfigurations,
        iupacLinkageBrackets,
        edgeNames,
        noEdgeLabels};
}

std::string cdsCondensedSequence::printSequence(const SequencePrintConfig& config, const SequenceData& sequence)
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
    SequencePrintData data = {
        config.derivativeBrackets,
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

std::string cdsCondensedSequence::printGraphViz(GraphVizDotConfig& configs, const SequenceData& sequence)
{
    std::vector<std::string> derivatives = sequenceDerivatives(sequence);
    std::vector<std::string> monosaccharideNames = sequenceMonosaccharideNames(sequence);
    graph::Graph graph = condensedSequenceGraph(sequence);

    std::vector<GraphVizResidueNode> nodes;
    nodes.reserve(nodeCount(graph));
    std::vector<GraphVizLinkage> linkages;
    linkages.reserve(edgeCount(graph));

    for (size_t n = 0; n < nodeCount(graph); n++)
    {
        size_t index = sourceNodeIndex(graph, n);
        std::string& monosaccharideName = monosaccharideNames[index];
        GraphVizImage image =
            findImage(configs, monosaccharideName, (sequence.residues.ringType[index] == "f") ? "f" : "");
        nodes.push_back({n, image, monosaccharideName, derivatives[index]});
    }

    for (size_t n = 0; n < edgeCount(graph); n++)
    {
        std::array<size_t, 2>& adj = graph.edges.nodeAdjacencies[n];
        size_t parent = adj[0];
        size_t child = adj[1];
        size_t childIndex = sourceNodeIndex(graph, child);
        std::string label = (configs.show_config_labels_ ? sequence.residues.configuration[childIndex] : "") +
                            (configs.show_position_labels_ ? edgeLinkage(sequence, graph.edges.indices[n]) : "");
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
        bool isAglycone = sequence.residues.type[nodeIndex] == cds::ResidueType::Aglycone;
        ss << (isAglycone ? graphVizAglyconeNode(nodes[n]) : graphVizSugarNode(nodes[n]));
    }
    for (auto& linkage : linkages)
    {
        ss << graphVizLinkageLine(nodes, linkage);
    }
    ss << "}\n";
    std::string str = ss.str();
    // Open and overwrite.
    codeUtils::writeToFile(configs.file_name_, [&str](std::ostream& stream) { stream << str; });
    return str;
}
