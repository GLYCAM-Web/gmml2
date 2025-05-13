#ifndef INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYGRAPH_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <array>
#include <vector>

namespace assembly
{
    struct Indices
    {
        size_t atomCount     = 0;
        size_t residueCount  = 0;
        size_t moleculeCount = 0;
        size_t assemblyCount = 0;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly;
    };

    struct Graph
    {
        Indices indices;
        graph::Graph atoms;
        graph::Graph residues;
        graph::Graph molecules;
        graph::Graph assemblies;
    };

    std::array<std::vector<bool>, 2> residueAtomsCloseToEdge(const assembly::Graph& graph, size_t residueEdgeIndex);
    std::vector<std::array<std::vector<bool>, 2>> atomsCloseToResidueEdges(const assembly::Graph& graph);

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
            return graph.indices.residueMolecule[graph.indices.atomResidue[n]] == moleculeId;
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

    inline size_t moleculeEdgeToResidueEdgeIndex(const Graph& graph, size_t edgeId)
    {
        return graph.molecules.edges.indices[edgeId];
    }

    inline size_t moleculeEdgeToAtomEdgeIndex(const Graph& graph, size_t edgeId)
    {
        return residueEdgeToAtomEdgeIndex(graph, moleculeEdgeToResidueEdgeIndex(graph, edgeId));
    }
} // namespace assembly

#endif
