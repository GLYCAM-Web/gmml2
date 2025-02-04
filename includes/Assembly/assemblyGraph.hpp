#ifndef INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP

#include "includes/Graph/graphTypes.hpp"

#include <vector>

namespace assembly
{
    struct Graph
    {
        size_t atomCount;
        size_t residueCount;
        size_t moleculeCount;
        size_t assemblyCount;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly;
        graph::Graph atoms;
        graph::Graph residues;
        graph::Graph molecules;
        graph::Graph assemblies;
    };

    inline const std::vector<size_t>& residueAtoms(const Graph& graph, size_t residueId)
    {
        return graph.residues.nodes.elements[residueId];
    };

    inline const std::vector<size_t>& moleculeResidues(const Graph& graph, size_t moleculeId)
    {
        return graph.molecules.nodes.elements[moleculeId];
    };

    inline const std::vector<size_t>& assemblyMolecules(const Graph& graph, size_t assemblyId)
    {
        return graph.assemblies.nodes.elements[assemblyId];
    };

    inline size_t residueEdgeToAtomEdgeIndex(const Graph& graph, size_t edgeId)
    {
        return graph.residues.edges.indices[edgeId];
    }

    inline size_t moleculeEdgeToAtomEdgeIndex(const Graph& graph, size_t edgeId)
    {
        return residueEdgeToAtomEdgeIndex(graph, graph.molecules.edges.indices[edgeId]);
    }
} // namespace assembly

#endif
