#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <string>
#include <vector>
#include <optional>
#include <variant>

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    auto residueLink = [&](size_t n)
    {
        return cdsCondensedSequence::getLink(sequence.nodes.type[n], sequence.nodes.linkage[n]);
    };
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        return residueLink(sequence.graph.edgeNodes[n][1]) > residueLink(sequence.graph.edgeNodes[k][1]);
    };

    return codeUtils::sortedBy(compare, edgeIds);
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::pruned(const SequenceData& sequence)
{
    size_t residueCount     = sequence.nodes.name.size();
    size_t initialEdgeCount = sequence.graph.edges.size();
    std::vector<size_t> aliveEdges;
    aliveEdges.reserve(initialEdgeCount);
    for (size_t n = 0; n < initialEdgeCount; n++)
    {
        const std::array<size_t, 2>& arr = sequence.graph.edgeNodes[n];
        if (sequence.graph.nodeAlive[arr[0]] && sequence.graph.nodeAlive[arr[1]])
        {
            aliveEdges.push_back(n);
        }
    }

    std::vector<size_t> residueOrder = codeUtils::boolsToIndices(sequence.graph.nodeAlive);

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

    for (size_t n : aliveEdges)
    {
        auto& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    NodeData nodes = {codeUtils::indicesToValues(sequence.nodes.fullString, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.type, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.name, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.linkage, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.chainPosition, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.ringType, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.configuration, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isomer, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.preIsomerModifier, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.ringShape, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.modifier, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isInternal, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isDerivative, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.nodeType, residueOrder)};

    EdgeData edges = {codeUtils::indicesToValues(sequence.edges.names, aliveEdges)};

    for (size_t n : residueOrder)
    {
        const NodeType& nodeType = sequence.nodes.nodeType[n];
        if (std::holds_alternative<ProbabilityNode>(nodeType) || std::holds_alternative<BranchNode>(nodeType))
        {
            throw std::runtime_error("Deactivate and bypass probability / branch nodes before manipulating sequence");
        }
    }

    return SequenceData {resultGraph, nodes, edges};
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::reordered(const SequenceData& sequence)
{
    size_t residueCount                               = sequence.nodes.name.size();
    std::vector<size_t> edgeIndexOrder                = codeUtils::indexVector(sequence.graph.edges);
    std::vector<size_t> edgeOrder                     = edgesSortedByLink(sequence, edgeIndexOrder);
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
        throw std::runtime_error("Sequence graph not fully connected");
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

    for (size_t n : edgeIndexOrder)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    NodeData nodes = {codeUtils::indicesToValues(sequence.nodes.fullString, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.type, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.name, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.linkage, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.chainPosition, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.ringType, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.configuration, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isomer, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.preIsomerModifier, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.ringShape, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.modifier, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isInternal, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.isDerivative, residueOrder),
                      codeUtils::indicesToValues(sequence.nodes.nodeType, residueOrder)};

    for (size_t n : residueOrder)
    {
        const NodeType& nodeType = sequence.nodes.nodeType[n];
        if (std::holds_alternative<ProbabilityNode>(nodeType) || std::holds_alternative<BranchNode>(nodeType))
        {
            throw std::runtime_error("Deactivate and bypass probability / branch nodes before manipulating sequence");
        }
    }

    return SequenceData {resultGraph, nodes, sequence.edges};
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::constructSequence(pcg32& rng, SequenceData sequence)
{
    graph::Graph traversable = graph::identity(sequence.graph);
    size_t nodeCount         = traversable.nodes.indices.size();
    std::vector<size_t> parent(nodeCount, 0);
    std::vector<size_t> closestLivingParent(nodeCount, 0);
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
        const cdsCondensedSequence::NodeType& nodeType = sequence.nodes.nodeType[current];
        // means we're dealing with a probability node
        if (std::holds_alternative<cdsCondensedSequence::ProbabilityNode>(nodeType))
        {
            const cdsCondensedSequence::ProbabilityNode& node =
                std::get<cdsCondensedSequence::ProbabilityNode>(nodeType);
            sequence.graph.nodeAlive[current] = false;
            included[current]   = codeUtils::uniformRandomDoubleWithinRange(rng, 0.0, 1.0) < node.probability;
            size_t choice       = codeUtils::weightedRandomIndex(rng, node.weights);
            size_t id           = node.heads[choice];
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
            const cdsCondensedSequence::BranchNode& node = std::get<cdsCondensedSequence::BranchNode>(nodeType);
            sequence.graph.nodeAlive[current]            = false;
            size_t id                                    = node.head;
            parent[id]                                   = current;
            traversed[id]                                = true;
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
    for (size_t n = 0; n < traversal.size(); n++)
    {
        size_t current = traversal[n];
        closestLivingParent[current] =
            included[parent[current]] ? parent[current] : closestLivingParent[parent[current]];
    }
    sequence.graph.nodeAlive = codeUtils::vectorAnd(sequence.graph.nodeAlive, included);

    // add missing connections
    for (size_t nodeId : codeUtils::boolsToIndices(included))
    {
        size_t n                                          = actualHead[nodeId];
        size_t closestParent                              = closestLivingParent[nodeId];
        cdsCondensedSequence::NodeType& nodeType          = sequence.nodes.nodeType[nodeId];
        cdsCondensedSequence::NodeType& closestParentType = sequence.nodes.nodeType[closestParent];
        if (std::holds_alternative<cdsCondensedSequence::BranchNode>(closestParentType))
        {
            if (std::holds_alternative<cdsCondensedSequence::BranchNode>(sequence.nodes.nodeType[nodeId]))
            {
                throw std::runtime_error("Head of branch dropped leaving nested branch with nowhere to connect: " +
                                         sequence.nodes.fullString[closestParent]);
            }
            isBranchHead[n]                          = true;
            cdsCondensedSequence::BranchNode& branch = std::get<cdsCondensedSequence::BranchNode>(closestParentType);
            std::string linkageName                  = sequence.nodes.linkage[n];
            size_t dashPosition                      = linkageName.find('-');
            if (dashPosition == std::string::npos)
            {
                throw std::runtime_error("Missing dash in branch linkage name: " + sequence.nodes.fullString[n]);
            }
            sequence.nodes.linkage[n] = linkageName.substr(0, dashPosition + 1) + branch.linkage;
            graph::addEdge(sequence.graph, {actualTail[closestLivingParent[closestParent]], n});
            sequence.edges.names.push_back("");
        }
        else if (std::holds_alternative<cdsCondensedSequence::ProbabilityNode>(nodeType) || !included[parent[nodeId]] ||
                 ((actualTail[closestParent] != parent[nodeId]) && (actualHead[closestParent] != n)))
        {
            graph::addEdge(sequence.graph, {actualTail[closestParent], n});
            sequence.edges.names.push_back("");
        }
    }

    for (size_t n = 0; n < sequence.graph.edges.size(); n++)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        size_t parent                     = edge[0];
        size_t child                      = edge[1];
        if (sequence.graph.nodeAlive[parent] && sequence.graph.nodeAlive[child])
        {
            bool isChildDeoxy = sequence.nodes.type[child] == cds::ResidueType::Deoxy;
            if (!isBranchHead[child] && !sequence.nodes.chainPosition[parent].empty())
            {
                std::string linkageName = sequence.nodes.linkage[child];
                size_t dashPosition     = linkageName.find('-');
                if (dashPosition != std::string::npos)
                {
                    sequence.nodes.linkage[child] =
                        linkageName.substr(0, dashPosition + 1) + sequence.nodes.chainPosition[parent];
                }
            }
            // It remains internal if it's already been made internal, or if the child is not a deoxy.
            sequence.nodes.isInternal[parent] = sequence.nodes.isInternal[parent] || !isChildDeoxy;
            sequence.edges.names[n]           = sequence.nodes.configuration[child] + sequence.nodes.linkage[child];
        }
    }
    return pruned(sequence);
}
