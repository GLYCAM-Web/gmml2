#ifndef INCLUDES_GRAPH_MANIPULATION_HPP
#define INCLUDES_GRAPH_MANIPULATION_HPP

#include "includes/Graph/types.hpp"

#include <cstddef>
#include <array>
#include <vector>

namespace graph
{
    size_t addNode(Database& graph);
    size_t addEdge(Database& graph, std::array<size_t, 2> nodes);
    void removeNode(Database& graph, size_t index);
    void removeEdge(Database& graph, size_t index);

    Database asData(const Graph& graph);
    Graph selectedQuotientAliveOrNot(const Database& graph, const std::vector<size_t>& nodeGroup,
                                     const std::vector<bool>& includeNode, const std::vector<bool>& includeEdge);
    Graph selectedQuotient(const Database& graph, const std::vector<size_t>& nodeGroup,
                           const std::vector<bool>& includeNode, const std::vector<bool>& includeEdge);
    Graph quotient(const Database& graph, const std::vector<size_t>& nodeGroup);
    Graph identity(const Database& graph);
} // namespace graph

#endif
