#ifndef INCLUDES_GRAPH_GRAPHFUNCTIONS_HPP
#define INCLUDES_GRAPH_GRAPHFUNCTIONS_HPP

#include "includes/Graph/graphTypes.hpp"

#include <cstddef>
#include <vector>

namespace graph
{
    std::vector<size_t> reachableNodes(const Graph& graph, const std::vector<bool>& excluded, size_t starting);
} // namespace graph

#endif
