#include "includes/Graph/graphManipulation.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <cstddef>
#include <array>
#include <vector>
#include <stdexcept>

namespace
{
    size_t groupMaxIndex(const std::vector<size_t>& group)
    {
        size_t maxIndex = 0;
        for (auto& n : group)
        {
            maxIndex = std::max(n, maxIndex);
        }
        return maxIndex;
    }

    std::pair<std::vector<size_t>, std::vector<size_t>> createNodes(const std::vector<size_t>& nodeGroup,
                                                                    const std::vector<bool>& includeNode)
    {
        std::vector<size_t> nodes;
        size_t groupMax = groupMaxIndex(nodeGroup);
        std::vector<size_t> groupIndex(groupMax + 1, -1);
        std::vector<bool> groupIndexSet(groupMax + 1, false);
        for (size_t n = 0; n < nodeGroup.size(); n++)
        {
            size_t group = nodeGroup[n];
            if (includeNode[n] && !groupIndexSet[group])
            {
                groupIndex[group]    = nodes.size();
                groupIndexSet[group] = true;
                nodes.push_back(group);
            }
        }
        return {nodes, groupIndex};
    }

    std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>>
    nodeAdjacencies(size_t nodeCount, const std::vector<std::array<size_t, 2>>& edgeNodes)
    {
        std::vector<std::vector<size_t>> nodes;
        nodes.resize(nodeCount);
        std::vector<std::vector<size_t>> edges;
        edges.resize(nodeCount);
        for (size_t n = 0; n < edgeNodes.size(); n++)
        {
            auto& nodePair = edgeNodes[n];
            size_t a       = nodePair[0];
            size_t b       = nodePair[1];
            nodes[a].push_back(b);
            nodes[b].push_back(a);
            edges[a].push_back(n);
            edges[b].push_back(n);
        }
        return {nodes, edges};
    }
} // namespace

size_t graph::addNode(Database& graph)
{
    size_t index = graph.nodes.size();
    graph.nodes.push_back(index);
    graph.nodeAlive.push_back(true);
    return index;
}

size_t graph::addEdge(Database& graph, std::array<size_t, 2> nodes)
{
    if (nodes[0] >= graph.nodes.size() || nodes[1] >= graph.nodes.size())
    {
        throw std::runtime_error("don't add connections between nodes that don't exist");
    }
    size_t index = graph.edges.size();
    graph.edges.push_back(index);
    graph.edgeAlive.push_back(true);
    graph.edgeNodes.push_back(nodes);
    return index;
}

void graph::removeNode(Database& graph, size_t index)
{
    if (index >= graph.nodeAlive.size())
    {
        throw std::runtime_error("don't remove nodes that don't exist");
    }
    graph.nodeAlive[index] = false;
}

void graph::removeEdge(Database& graph, size_t index)
{
    if (index >= graph.edgeAlive.size())
    {
        throw std::runtime_error("don't remove edges that don't exist");
    }
    graph.edgeAlive[index] = false;
}

void graph::reserveNodes(Database& graph, size_t size)
{
    graph.nodes.reserve(size);
    graph.nodeAlive.reserve(size);
}

void graph::reserveEdges(Database& graph, size_t size)
{
    graph.edges.reserve(size);
    graph.edgeNodes.reserve(size);
    graph.edgeAlive.reserve(size);
}

graph::Database graph::asData(const Graph& graph)
{
    std::vector<size_t> nodes = codeUtils::indexVector(graph.nodes.indices);
    std::vector<bool> nodeAlive(nodes.size(), true);
    std::vector<size_t> edges = codeUtils::indexVector(graph.edges.indices);
    std::vector<bool> edgeAlive(edges.size(), true);
    return {nodes, nodeAlive, edges, edgeAlive, graph.edges.nodeAdjacencies};
}

graph::Graph graph::selectedQuotientAliveOrNot(const Database& graph, const std::vector<size_t>& nodeGroup,
                                               const std::vector<bool>& includeNode,
                                               const std::vector<bool>& includeEdge)
{
    auto nodeData                  = createNodes(nodeGroup, includeNode);
    std::vector<size_t> nodes      = nodeData.first;
    std::vector<size_t> groupIndex = nodeData.second;
    std::vector<std::vector<size_t>> elements(nodes.size(), std::vector<size_t> {});
    for (size_t n = 0; n < graph.nodes.size(); n++)
    {
        if (includeNode[n])
        {
            elements[groupIndex[nodeGroup[n]]].push_back(graph.nodes[n]);
        }
    }
    std::vector<size_t> edges;
    std::vector<std::array<size_t, 2>> edgeNodes;
    for (size_t n = 0; n < graph.edges.size(); n++)
    {
        auto& en  = graph.edgeNodes[n];
        size_t a  = en[0];
        size_t b  = en[1];
        size_t ag = nodeGroup[a];
        size_t bg = nodeGroup[b];
        if (includeEdge[n] && includeNode[a] && includeNode[b] && (ag != bg))
        {
            edges.push_back(graph.edges[n]);
            edgeNodes.push_back({groupIndex[ag], groupIndex[bg]});
        }
    }
    auto adjacencies                           = nodeAdjacencies(nodes.size(), edgeNodes);
    std::vector<std::vector<size_t>> nodeNodes = adjacencies.first;
    std::vector<std::vector<size_t>> nodeEdges = adjacencies.second;
    return {
        {nodes, nodeNodes, nodeEdges, elements},
        {edges, edgeNodes},
        graph
    };
}

graph::Graph graph::selectedQuotient(const Database& graph, const std::vector<size_t>& nodeGroup,
                                     const std::vector<bool>& includeNode, const std::vector<bool>& includeEdge)
{
    return selectedQuotientAliveOrNot(graph, nodeGroup, codeUtils::vectorAnd(graph.nodeAlive, includeNode),
                                      codeUtils::vectorAnd(graph.edgeAlive, includeEdge));
}

graph::Graph graph::subgraph(const Database& graph, const std::vector<bool>& includeNode,
                             const std::vector<bool>& includeEdge)
{
    return selectedQuotient(graph, graph.nodes, includeNode, includeEdge);
}

graph::Graph graph::quotient(const Database& graph, const std::vector<size_t>& nodeGroup)
{
    return selectedQuotient(graph, nodeGroup, graph.nodeAlive, graph.edgeAlive);
}

graph::Graph graph::identity(const Database& graph)
{
    return quotient(graph, graph.nodes);
}
