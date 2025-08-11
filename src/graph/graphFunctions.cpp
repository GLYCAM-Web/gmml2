#include "include/graph/graphFunctions.hpp"

#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"

#include <cstddef>
#include <vector>

namespace gmml
{
    namespace graph
    {
        bool edgeAlive(const Database& db, size_t edgeId)
        {
            const std::array<size_t, 2>& nodes = db.edges.nodes[edgeId];
            return db.nodes.alive[nodes[0]] && db.nodes.alive[nodes[1]];
        }

        std::vector<bool> reachableNodes(const Graph& graph, const std::vector<bool>& excluded, size_t starting)
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
                    if (!(excluded[adj] || reachable[adj] || util::contains(toVisit, adj)))
                    {
                        reachable[adj] = true;
                        toVisit.push_back(adj);
                    }
                }
            }
            return reachable;
        }
    } // namespace graph
} // namespace gmml
