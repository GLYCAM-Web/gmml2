#include "includes/Graph/graphFunctions.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <cstddef>
#include <vector>

std::vector<size_t> graph::reachableNodes(const Graph& graph, const std::vector<bool>& excluded, size_t starting)
{
    std::vector<bool> reachable(nodeCount(graph), false);
    std::vector<size_t> toVisit;
    toVisit.reserve(nodeCount(graph));
    toVisit.push_back(starting);
    reachable[starting] = !excluded[starting];
    while (!toVisit.empty())
    {
        size_t visiting = toVisit.back();
        toVisit.pop_back();
        for (size_t adj : graph.nodes.nodeAdjacencies[visiting])
        {
            if (!(excluded[adj] || reachable[adj] || codeUtils::contains(toVisit, adj)))
            {
                reachable[adj] = true;
                toVisit.push_back(adj);
            }
        }
    }
    return codeUtils::boolsToIndices(reachable);
}
