#ifndef INCLUDES_GRAPH_TYPES_HPP
#define INCLUDES_GRAPH_TYPES_HPP

#include <cstddef>
#include <array>
#include <vector>

namespace graph
{
    struct Database
    {
        Database(const std::vector<size_t>& nodes_, const std::vector<bool>& nodeAlive_,
                 const std::vector<size_t>& edges_, const std::vector<bool>& edgeAlive_,
                 const std::vector<std::array<size_t, 2>>& edgeNodes_)
            : nodes(nodes_), nodeAlive(nodeAlive_), edges(edges_), edgeAlive(edgeAlive_), edgeNodes(edgeNodes_)
        {}

        std::vector<size_t> nodes;
        std::vector<bool> nodeAlive;
        std::vector<size_t> edges;
        std::vector<bool> edgeAlive;
        std::vector<std::array<size_t, 2>> edgeNodes;
    };

    struct GraphNodes
    {
        GraphNodes(const std::vector<size_t>& indices_, const std::vector<std::vector<size_t>>& nodeAdjacencies_,
                   const std::vector<std::vector<size_t>>& edgeAdjacencies_,
                   const std::vector<std::vector<size_t>>& elements_)
            : indices(indices_), nodeAdjacencies(nodeAdjacencies_), edgeAdjacencies(edgeAdjacencies_),
              elements(elements_)
        {}

        std::vector<size_t> indices;
        std::vector<std::vector<size_t>> nodeAdjacencies;
        std::vector<std::vector<size_t>> edgeAdjacencies;
        std::vector<std::vector<size_t>> elements;
    };

    struct GraphEdges
    {
        GraphEdges(const std::vector<size_t>& indices_, const std::vector<std::array<size_t, 2>>& nodeAdjacencies_)
            : indices(indices_), nodeAdjacencies(nodeAdjacencies_)
        {}

        std::vector<size_t> indices;
        std::vector<std::array<size_t, 2>> nodeAdjacencies;
    };

    struct Graph
    {
        Graph(const GraphNodes& nodes_, const GraphEdges& edges_, const Database& source_)
            : nodes(nodes_), edges(edges_), source(source_)
        {}

        GraphNodes nodes;
        GraphEdges edges;
        Database source;
    };
} // namespace graph

#endif
