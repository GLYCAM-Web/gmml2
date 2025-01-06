#ifndef INCLUDES_GRAPH_GRAPHMATCHING_HPP
#define INCLUDES_GRAPH_GRAPHMATCHING_HPP

#include "includes/Graph/graphTypes.hpp"

#include <cstddef>
#include <string>
#include <array>
#include <vector>

namespace graph
{
    struct LabeledGraph
    {
        const Database& db;
        const std::vector<std::string>& names;
    };

    bool graphsMatch(const LabeledGraph& reference, const LabeledGraph& other);
} // namespace graph

#endif
