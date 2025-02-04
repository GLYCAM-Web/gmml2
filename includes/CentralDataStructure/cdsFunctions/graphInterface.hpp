#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <vector>

namespace cds
{
    struct GraphIndexData
    {
        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules;
        std::vector<Assembly*> assemblies;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly;
    };

    GraphIndexData toIndexData(const std::vector<Residue*> inputResidues);
    GraphIndexData toIndexData(const std::vector<Molecule*> molecules);
    GraphIndexData toIndexData(const std::vector<Assembly*> assemblies);
    graph::Database createGraphData(const GraphIndexData& indices);
    assembly::Graph createAssemblyGraph(const GraphIndexData& indices, const std::vector<bool>& includedAtoms);
} // namespace cds

#endif
