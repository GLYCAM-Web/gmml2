#ifndef INCLUDE_GRAPH_GRAPHTYPES_HPP
#define INCLUDE_GRAPH_GRAPHTYPES_HPP

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
            std::vector<size_t> localIndices;
            std::vector<size_t> sourceIndices;
            std::vector<std::vector<size_t>> nodeAdjacencies;
            std::vector<std::vector<size_t>> edgeAdjacencies;
            std::vector<std::vector<size_t>> constituents;
        };

        struct GraphEdges
        {
            std::vector<size_t> localIndices;
            std::vector<size_t> sourceIndices;
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

        inline size_t nodeCount(const Graph& graph) { return graph.nodes.localIndices.size(); }

        inline size_t edgeCount(const Graph& graph) { return graph.edges.localIndices.size(); }

        inline const std::vector<size_t>& localIndices(const GraphNodes& nodes) { return nodes.localIndices; }

        inline const std::vector<size_t>& sourceIndices(const GraphNodes& nodes) { return nodes.sourceIndices; }

        inline const std::vector<size_t>& localIndices(const GraphEdges& edges) { return edges.localIndices; }

        inline const std::vector<size_t>& sourceIndices(const GraphEdges& edges) { return edges.sourceIndices; }

        inline size_t sourceNodeIndex(const Graph& graph, size_t n) { return graph.nodes.sourceIndices[n]; }

        inline size_t sourceEdgeIndex(const Graph& graph, size_t n) { return graph.edges.sourceIndices[n]; }
    } // namespace graph
} // namespace gmml

#endif
