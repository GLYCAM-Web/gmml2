#ifndef INCLUDES_GRAPH_GRAPHTYPES_HPP
#define INCLUDES_GRAPH_GRAPHTYPES_HPP

#include <cstddef>
#include <array>
#include <vector>

namespace graph
{
    struct Database
    {
        std::vector<size_t> nodes;
        std::vector<bool> nodeAlive;
        std::vector<size_t> edges;
        std::vector<bool> edgeAlive;
        std::vector<std::array<size_t, 2>> edgeNodes;
    };

    struct GraphNodes
    {
        std::vector<size_t> indices;
        std::vector<std::vector<size_t>> nodeAdjacencies;
        std::vector<std::vector<size_t>> edgeAdjacencies;
        std::vector<std::vector<size_t>> elements;
    };

    struct GraphEdges
    {
        std::vector<size_t> indices;
        std::vector<std::array<size_t, 2>> nodeAdjacencies;
    };

    struct Graph
    {
        GraphNodes nodes;
        GraphEdges edges;
        Database source;
    };
} // namespace graph

#endif
