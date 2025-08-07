#include "include/assembly/assemblyGraph.hpp"

#include "include/graph/graphManipulation.hpp"
#include "include/util/containers.hpp"

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

        std::vector<size_t> moleculeAtoms(const Graph& graph, size_t moleculeId)
        {
            std::function<bool(const size_t&)> inMolecule = [&](const size_t& n)
            { return graph.indices.residueMolecule[graph.indices.atomResidue[n]] == moleculeId; };
            std::vector<bool> selected = util::vectorMap(inMolecule, sourceIndices(graph.atoms.nodes));
            return util::boolsToIndices(selected);
        }

        Graph createAssemblyGraph(const Indices& indices, const graph::Database& atomGraphData)
        {
            graph::Graph atomGraph = graph::identity(atomGraphData);
            graph::Graph residueGraph = graph::quotient(atomGraphData, indices.atomResidue);
            graph::Database residueData = graph::asData(residueGraph);
            graph::Graph moleculeGraph =
                graph::quotient(residueData, util::indicesToValues(indices.residueMolecule, residueData.nodes));
            graph::Database moleculeData = graph::asData(moleculeGraph);
            graph::Graph assemblyGraph =
                graph::quotient(moleculeData, util::indicesToValues(indices.moleculeAssembly, moleculeData.nodes));
            return assembly::Graph {
                {indices.atomCount,
                 indices.residueCount,
                 indices.moleculeCount,
                 indices.assemblyCount,
                 atomGraphData.nodeAlive,
                 indices.atomResidue,
                 indices.residueMolecule,
                 indices.moleculeAssembly},
                atomGraph,
                residueGraph,
                moleculeGraph,
                assemblyGraph
            };
        }
    } // namespace assembly
} // namespace gmml
