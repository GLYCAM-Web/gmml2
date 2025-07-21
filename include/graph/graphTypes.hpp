#ifndef INCLUDES_GRAPH_GRAPHTYPES_HPP
#define INCLUDES_GRAPH_GRAPHTYPES_HPP

#include <array>
#include <cstddef>
#include <vector>

namespace gmml
{
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
            std::vector<std::vector<size_t>> constituents;
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

        inline size_t nodeCount(const Database& db) { return db.nodes.size(); }

        inline size_t edgeCount(const Database& db) { return db.edges.size(); }

        inline size_t nodeCount(const Graph& graph) { return graph.nodes.indices.size(); }

        inline size_t edgeCount(const Graph& graph) { return graph.edges.indices.size(); }

        inline size_t sourceNodeIndex(const Graph& graph, size_t n) { return graph.nodes.indices[n]; }
    } // namespace graph
} // namespace gmml

#endif
