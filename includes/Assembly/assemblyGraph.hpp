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
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        graph::Graph atoms;
        graph::Graph residues;
        graph::Graph molecules;
    };

    inline const std::vector<size_t>& residueAtoms(const Graph& graph, size_t residueId)
    {
        return graph.residues.nodes.elements[residueId];
    };

    inline const std::vector<size_t>& moleculeResidues(const Graph& graph, size_t moleculeId)
    {
        return graph.molecules.nodes.elements[moleculeId];
    };

    inline size_t moleculeEdgeToAtomEdgeIndex(const Graph& graph, size_t moleculeEdgeId)
    {
        size_t residueBondIndex = graph.molecules.edges.indices[moleculeEdgeId];
        return graph.residues.edges.indices[residueBondIndex];
    }
} // namespace assembly

#endif
