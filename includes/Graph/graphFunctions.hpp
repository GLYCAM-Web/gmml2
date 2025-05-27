#ifndef INCLUDES_GRAPH_GRAPHFUNCTIONS_HPP
#define INCLUDES_GRAPH_GRAPHFUNCTIONS_HPP

#include "includes/Graph/graphTypes.hpp"

#include <cstddef>
#include <vector>

namespace graph
{
    bool edgeAlive(const Database& db, size_t edgeId);
    std::vector<bool> reachableNodes(const Graph& graph, const std::vector<bool>& excluded, size_t starting);
} // namespace graph

#endif
