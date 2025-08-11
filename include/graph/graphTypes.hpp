#ifndef INCLUDE_GRAPH_GRAPHTYPES_HPP
#define INCLUDE_GRAPH_GRAPHTYPES_HPP

#include <array>
#include <cstddef>
#include <vector>

namespace gmml
{
    namespace graph
    {
        struct NodeDB
        {
            std::vector<size_t> indices;
            std::vector<bool> alive;
        };

        struct EdgeDB
        {
            std::vector<size_t> indices;
            std::vector<bool> alive;
            std::vector<std::array<size_t, 2>> nodes;
        };

        struct Database
        {
            NodeDB nodes;
            EdgeDB edges;
        };

        struct GraphNodes
        {
            std::vector<size_t> aliveIndices;
            std::vector<std::vector<size_t>> nodeAdjacencies;
            std::vector<std::vector<size_t>> edgeAdjacencies;
            std::vector<std::vector<size_t>> constituents;
        };

        struct GraphEdges
        {
            std::vector<size_t> sourceIndices;
            std::vector<std::array<size_t, 2>> nodeAdjacencies;
        };

        struct Graph
        {
            GraphNodes nodes;
            GraphEdges edges;
            Database source;
        };

        inline size_t nodeCount(const Database& db) { return db.nodes.indices.size(); }

        inline size_t edgeCount(const Database& db) { return db.edges.indices.size(); }

        inline size_t nodeCount(const Graph& graph) { return graph.nodes.aliveIndices.size(); }

        inline size_t edgeCount(const Graph& graph) { return graph.edges.sourceIndices.size(); }

        inline size_t sourceNodeCount(const Graph& graph) { return graph.nodes.constituents.size(); }

        inline const std::vector<size_t>& indices(const GraphNodes& nodes) { return nodes.aliveIndices; }

        inline size_t sourceIndex(const GraphEdges& edges, size_t n) { return edges.sourceIndices[n]; }

    } // namespace graph
} // namespace gmml

#endif
