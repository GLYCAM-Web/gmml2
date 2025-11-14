#include "include/graph/graphManipulation.hpp"

#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace gmml
{
    namespace graph
    {
        namespace
        {
            std::vector<bool> aliveGroups(
                size_t groupCount, const std::vector<size_t>& nodeGroup, const std::vector<bool>& includeNode)
            {
                std::vector<bool> groupAlive(groupCount, false);
                for (size_t group : util::boolsToValues(nodeGroup, includeNode))
                {
                    groupAlive[group] = true;
                }
                return groupAlive;
            }

            std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> nodeAdjacencies(
                size_t nodeCount, const std::vector<std::array<size_t, 2>>& edgeNodes)
            {
                std::vector<std::vector<size_t>> nodes;
                nodes.resize(nodeCount);
                std::vector<std::vector<size_t>> edges;
                edges.resize(nodeCount);
                for (size_t n = 0; n < edgeNodes.size(); n++)
                {
                    auto& nodePair = edgeNodes[n];
                    size_t a = nodePair[0];
                    size_t b = nodePair[1];
                    nodes[a].push_back(b);
                    nodes[b].push_back(a);
                    edges[a].push_back(n);
                    edges[b].push_back(n);
                }
                return {nodes, edges};
            }
        } // namespace

        size_t addNode(Database& graph)
        {
            size_t index = graph.nodes.indices.size();
            graph.nodes.indices.push_back(index);
            graph.nodes.alive.push_back(true);
            return index;
        }

        size_t addEdge(Database& graph, std::array<size_t, 2> nodes)
        {
            if (nodes[0] >= graph.nodes.indices.size() || nodes[1] >= graph.nodes.indices.size())
            {
                throw std::runtime_error("don't add connections between nodes that don't exist");
            }
            size_t index = graph.edges.indices.size();
            graph.edges.indices.push_back(index);
            graph.edges.alive.push_back(true);
            graph.edges.nodes.push_back(nodes);
            return index;
        }

        void removeNode(Database& graph, size_t index)
        {
            if (index >= graph.nodes.alive.size())
            {
                throw std::runtime_error("don't remove nodes that don't exist");
            }
            graph.nodes.alive[index] = false;
        }

        void removeEdge(Database& graph, size_t index)
        {
            if (index >= graph.edges.alive.size())
            {
                throw std::runtime_error("don't remove edges that don't exist");
            }
            graph.edges.alive[index] = false;
        }

        void reserveNodes(Database& graph, size_t size)
        {
            graph.nodes.indices.reserve(size);
            graph.nodes.alive.reserve(size);
        }

        void reserveEdges(Database& graph, size_t size)
        {
            graph.edges.indices.reserve(size);
            graph.edges.nodes.reserve(size);
            graph.edges.alive.reserve(size);
        }

        Database asData(const Graph& graph)
        {
            std::vector<bool> nodeAlive =
                util::indicesToBools(graph.nodes.constituents.size(), graph.nodes.aliveIndices);
            std::vector<bool> edgeAlive(edgeCount(graph), true);
            return {
                {util::indexVector(nodeAlive), nodeAlive},
                {util::indexVector(edgeAlive), edgeAlive, graph.edges.nodeAdjacencies}
            };
        }

        Graph selectedQuotientAliveOrNot(
            const Database& source,
            size_t groupCount,
            const std::vector<size_t>& nodeGroup,
            const std::vector<bool>& includeNode,
            const std::vector<bool>& includeEdge)
        {
            std::vector<bool> groupAlive = aliveGroups(groupCount, nodeGroup, includeNode);
            std::vector<std::vector<size_t>> constituents(groupCount, std::vector<size_t> {});
            for (size_t n = 0; n < source.nodes.indices.size(); n++)
            {
                if (includeNode[n])
                {
                    constituents[nodeGroup[n]].push_back(source.nodes.indices[n]);
                }
            }
            std::vector<size_t> edges;
            std::vector<std::array<size_t, 2>> edgeNodes;
            for (size_t n : source.edges.indices)
            {
                const std::array<size_t, 2>& en = source.edges.nodes[n];
                size_t a = en[0];
                size_t b = en[1];
                size_t ag = nodeGroup[a];
                size_t bg = nodeGroup[b];
                if (includeEdge[n] && includeNode[a] && includeNode[b] && (ag != bg))
                {
                    edges.push_back(source.edges.indices[n]);
                    edgeNodes.push_back({ag, bg});
                }
            }
            auto adjacencies = nodeAdjacencies(groupCount, edgeNodes);
            std::vector<std::vector<size_t>>& nodeNodes = adjacencies.first;
            std::vector<std::vector<size_t>>& nodeEdges = adjacencies.second;
            return {
                {util::indexVector(groupCount), util::boolsToIndices(groupAlive), nodeNodes, nodeEdges, constituents},
                {edges, edgeNodes},
                source
            };
        }

        Graph selectedQuotient(
            const Database& source,
            size_t groupCount,
            const std::vector<size_t>& nodeGroup,
            const std::vector<bool>& includeNode,
            const std::vector<bool>& includeEdge)
        {
            return selectedQuotientAliveOrNot(
                source,
                groupCount,
                nodeGroup,
                util::vectorAnd(source.nodes.alive, includeNode),
                util::vectorAnd(source.edges.alive, includeEdge));
        }

        Graph subgraph(
            const Database& source, const std::vector<bool>& includeNode, const std::vector<bool>& includeEdge)
        {
            return selectedQuotient(
                source, source.nodes.indices.size(), source.nodes.indices, includeNode, includeEdge);
        }

        Graph quotient(const Database& source, size_t groupCount, const std::vector<size_t>& nodeGroup)
        {
            return selectedQuotient(source, groupCount, nodeGroup, source.nodes.alive, source.edges.alive);
        }

        Graph identity(const Database& source)
        {
            return quotient(source, source.nodes.indices.size(), source.nodes.indices);
        }

    } // namespace graph
} // namespace gmml
