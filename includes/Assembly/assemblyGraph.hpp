#ifndef INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/CodeUtils/containers.hpp"

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
        return graph.residues.nodes.constituents[residueId];
    };

    inline const std::vector<size_t>& moleculeResidues(const Graph& graph, size_t moleculeId)
    {
        return graph.molecules.nodes.constituents[moleculeId];
    };

    inline std::vector<size_t> moleculeAtoms(const Graph& graph, size_t moleculeId)
    {
        std::function<bool(const size_t&)> inMolecule = [&](const size_t& n)
        {
            return graph.residueMolecule[graph.atomResidue[n]] == moleculeId;
        };
        std::vector<bool> selected = codeUtils::vectorMap(inMolecule, graph.atoms.nodes.indices);
        return codeUtils::boolsToIndices(selected);
    }

    inline const std::vector<size_t>& assemblyMolecules(const Graph& graph, size_t assemblyId)
    {
        return graph.assemblies.nodes.constituents[assemblyId];
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
