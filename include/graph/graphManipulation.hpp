#ifndef INCLUDE_GRAPH_GRAPHMANIPULATION_HPP
#define INCLUDE_GRAPH_GRAPHMANIPULATION_HPP

#include "include/graph/graphTypes.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace gmml
{
    namespace graph
    {
        size_t addNode(Database& graph);
        size_t addEdge(Database& graph, std::array<size_t, 2> nodes);
        void removeNode(Database& graph, size_t index);
        void removeEdge(Database& graph, size_t index);

        void reserveNodes(Database& graph, size_t size);
        void reserveEdges(Database& graph, size_t size);

        Database asData(const Graph& graph);
        Graph selectedQuotientAliveOrNot(
            const Database& graph,
            const std::vector<size_t>& nodeGroup,
            const std::vector<bool>& includeNode,
            const std::vector<bool>& includeEdge);
        Graph selectedQuotient(
            const Database& graph,
            const std::vector<size_t>& nodeGroup,
            const std::vector<bool>& includeNode,
            const std::vector<bool>& includeEdge);
        Graph subgraph(
            const Database& graph, const std::vector<bool>& includeNode, const std::vector<bool>& includeEdge);
        Graph quotient(const Database& graph, const std::vector<size_t>& nodeGroup);
        Graph identity(const Database& graph);
    } // namespace graph
} // namespace gmml

#endif
