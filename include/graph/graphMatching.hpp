#ifndef INCLUDES_GRAPH_GRAPHMATCHING_HPP
#define INCLUDES_GRAPH_GRAPHMATCHING_HPP

#include "include/graph/graphTypes.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace gmml
{
    namespace graph
    {
        struct LabeledGraph
        {
            const Database& db;
            const std::vector<std::string>& names;
        };

        bool graphsMatch(const LabeledGraph& reference, const LabeledGraph& other);
    } // namespace graph
} // namespace gmml

#endif
