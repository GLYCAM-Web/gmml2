#include "include/assembly/assemblyGraph.hpp"

#include <array>
#include <vector>

namespace gmml
{
    namespace
    {
        std::vector<bool> atomsCloseToBond(const assembly::Graph& graph, size_t residueIndex, size_t atomIndex)
        {
            const std::vector<size_t>& atomAdj = graph.atoms.nodes.nodeAdjacencies[atomIndex];
            const std::vector<size_t>& atomIndices = residueAtoms(graph, residueIndex);
            std::vector<bool> close(atomIndices.size(), false);
            close[util::indexOf(atomIndices, atomIndex)] = true;
            for (size_t adj : atomAdj)
            {
                if (graph.indices.atomResidue[adj] == residueIndex)
                {
                    close[util::indexOf(atomIndices, adj)] = true;
                }
            }
            return close;
        };
    } // namespace

    namespace assembly
    {
        std::array<std::vector<bool>, 2> residueAtomsCloseToEdge(const assembly::Graph& graph, size_t residueEdgeIndex)
        {
            size_t atomEdgeIndex = residueEdgeToAtomEdgeIndex(graph, residueEdgeIndex);
            const std::array<size_t, 2>& residues = graph.residues.edges.nodeAdjacencies[residueEdgeIndex];
            const std::array<size_t, 2>& atoms = graph.atoms.edges.nodeAdjacencies[atomEdgeIndex];
            return {atomsCloseToBond(graph, residues[0], atoms[0]), atomsCloseToBond(graph, residues[1], atoms[1])};
        }

        std::vector<std::array<std::vector<bool>, 2>> atomsCloseToResidueEdges(const assembly::Graph& graph)
        {
            size_t edgeCount = graph::edgeCount(graph.residues);
            std::vector<std::array<std::vector<bool>, 2>> result;
            result.reserve(edgeCount);
            for (size_t n = 0; n < edgeCount; n++)
            {
                result.push_back(residueAtomsCloseToEdge(graph, n));
            }
            return result;
        }
    } // namespace assembly
} // namespace gmml
