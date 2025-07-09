#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <stdexcept>

std::string cdsCondensedSequence::mainLinkage(cds::ResidueType type, const std::string& linkage)
{
    switch (type)
    {
        case (cds::ResidueType::Sugar):
            return linkage.substr(linkage.size() - 1, 1);
        case (cds::ResidueType::Derivative):
            return linkage.substr(0, 1);
        case (cds::ResidueType::Deoxy):
            return linkage.substr(0, 1);
        default:
            return std::string("0");
    }
}

std::string cdsCondensedSequence::mainLinkage(const SequenceData& sequence, size_t residueId)
{
    return mainLinkage(sequence.residues.type[residueId], sequence.residues.linkage[residueId]);
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::instantiate(const AbstractSequence& data)
{
    SequenceData result;
    for (size_t n = 0; n < data.nodes.size(); n++)
    {
        const ResidueNode& node                   = data.nodes[n];
        const ParsedResidueComponents& components = node.components;
        graph::addNode(result.graph);
        result.residues.fullString.push_back(components.fullString);
        result.residues.type.push_back(components.type);
        result.residues.name.push_back(components.name);
        result.residues.linkage.push_back(components.linkage);
        result.residues.ringType.push_back(components.ringType);
        result.residues.configuration.push_back(components.configuration);
        result.residues.isomer.push_back(components.isomer);
        result.residues.preIsomerModifier.push_back(components.preIsomerModifier);
        result.residues.ringShape.push_back(components.ringShape);
        result.residues.modifier.push_back(components.modifier);
        result.residues.isInternal.push_back(false);
        result.residues.isDerivative.push_back(
            codeUtils::contains({cds::ResidueType::Deoxy, cds::ResidueType::Derivative}, components.type));
    }
    for (auto& edge : data.graph.edgeNodes)
    {
        size_t parent = edge[0];
        size_t child  = edge[1];

        bool isSugar         = result.residues.type[child] == cds::ResidueType::Sugar;
        std::string edgeName = (isSugar ? result.residues.configuration[child] : "") + result.residues.linkage[child];
        graph::addEdge(result.graph, edge);
        result.edges.names.push_back(edgeName);
        // It remains internal if it's already been made internal, or if the child is not a deoxy.
        bool isChildDeoxy                  = result.residues.type[child] == cds::ResidueType::Deoxy;
        result.residues.isInternal[parent] = (!isChildDeoxy || result.residues.isInternal[parent]);
    }
    return result;
}

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        if (sequence.graph.edgeNodes[n][0] == sequence.graph.edgeNodes[k][0])
        {
            return mainLinkage(sequence, sequence.graph.edgeNodes[n][1]) >
                   mainLinkage(sequence, sequence.graph.edgeNodes[k][1]);
        }
        else
        {
            return sequence.graph.edgeNodes[n][0] < sequence.graph.edgeNodes[k][0];
        }
    };

    return codeUtils::sortedBy(compare, edgeIds);
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

    for (size_t n : edgeOrder)
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
                            codeUtils::indicesToValues(sequence.residues.preIsomerModifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringShape, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.modifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isInternal, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isDerivative, residueOrder)};

    EdgeData edges = {codeUtils::indicesToValues(sequence.edges.names, edgeOrder)};

    return SequenceData {resultGraph, residues, edges};
}
