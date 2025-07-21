#include "include/graph/graphMatching.hpp"

#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace gmml
{
    namespace graph
    {
        namespace
        {
            bool equalSize(const Database& a, const Database& b)
            {
                return (a.nodes.size() == b.nodes.size()) && (a.edges.size() == b.edges.size());
            }
        } // namespace

        bool graphsMatch(const LabeledGraph& reference, const LabeledGraph& other)
        {
            if (!equalSize(reference.db, other.db))
            {
                return false;
            }
            else
            {
                // we assign indices of the reference nodes to the other nodes
                // if they can be assigned a unique index, the nodes match
                std::vector<size_t> nodeIndices = util::indexVector(reference.names);
                std::vector<size_t> nodeOrder;
                nodeOrder.reserve(reference.names.size());
                for (auto& name : other.names)
                {
                    auto nameEquals = [&](size_t n) { return reference.names[n] == name; };
                    auto it = std::find_if(nodeIndices.begin(), nodeIndices.end(), nameEquals);
                    if (it == nodeIndices.end())
                    {
                        return false;
                    }
                    nodeOrder.push_back(*it);
                    nodeIndices.erase(it);
                }
                // same for the edges, they match if they can be assigned an index
                // unlike nodes, we don't need to store the assigned indices explicitly
                std::vector<size_t> edgeIndices = util::indexVector(reference.db.edges);
                for (auto& adj : other.db.edgeNodes)
                {
                    auto edgeEquals = [&](size_t n)
                    {
                        const std::array<size_t, 2>& nodes = reference.db.edgeNodes[n];
                        auto dirEqual = [&](int dir)
                        { return nodeOrder[adj[0]] == nodes[0 != dir] && nodeOrder[adj[1]] == nodes[1 != dir]; };
                        return dirEqual(0) || dirEqual(1);
                    };
                    auto it = std::find_if(edgeIndices.begin(), edgeIndices.end(), edgeEquals);
                    if (it == edgeIndices.end())
                    {
                        return false;
                    }
                    edgeIndices.erase(it);
                }
                return true;
            }
        }
    } // namespace graph
} // namespace gmml
