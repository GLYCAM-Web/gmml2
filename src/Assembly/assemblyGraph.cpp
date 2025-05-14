#include "includes/Assembly/assemblyGraph.hpp"

#include <array>
#include <vector>

namespace
{
    std::vector<bool> atomsCloseToBond(const assembly::Graph& graph, size_t residueIndex, size_t atomIndex)
    {
        const std::vector<size_t>& atomAdj     = graph.atoms.nodes.nodeAdjacencies[atomIndex];
        const std::vector<size_t>& atomIndices = residueAtoms(graph, residueIndex);
        std::vector<bool> close(atomIndices.size(), false);
        close[codeUtils::indexOf(atomIndices, atomIndex)] = true;
        for (size_t adj : atomAdj)
        {
            if (graph.indices.atomResidue[adj] == residueIndex)
            {
                close[codeUtils::indexOf(atomIndices, adj)] = true;
            }
        }
        return close;
    };
} // namespace

std::array<std::vector<bool>, 2> assembly::residueAtomsCloseToEdge(const assembly::Graph& graph,
                                                                   size_t residueEdgeIndex)
{
    size_t atomEdgeIndex                  = residueEdgeToAtomEdgeIndex(graph, residueEdgeIndex);
    const std::array<size_t, 2>& residues = graph.residues.edges.nodeAdjacencies[residueEdgeIndex];
    const std::array<size_t, 2>& atoms    = graph.atoms.edges.nodeAdjacencies[atomEdgeIndex];
    return {atomsCloseToBond(graph, residues[0], atoms[0]), atomsCloseToBond(graph, residues[1], atoms[1])};
}

std::vector<std::array<std::vector<bool>, 2>> assembly::atomsCloseToResidueEdges(const assembly::Graph& graph)
{
    size_t edgeCount = graph.residues.edges.indices.size();
    std::vector<std::array<std::vector<bool>, 2>> result;
    result.reserve(edgeCount);
    for (size_t n = 0; n < edgeCount; n++)
    {
        result.push_back(residueAtomsCloseToEdge(graph, n));
    }
    return result;
}
