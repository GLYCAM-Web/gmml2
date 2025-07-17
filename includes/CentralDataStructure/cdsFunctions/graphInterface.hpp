#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/CodeUtils/containerTypes.hpp"

#include <vector>

namespace cds
{
    struct GraphObjects
    {
        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules;
        std::vector<Assembly*> assemblies;
    };

    struct GraphIndexData
    {
        assembly::Indices indices;
        GraphObjects objects;
    };

    GraphIndexData toIndexData(const std::vector<Residue*> inputResidues);
    GraphIndexData toIndexData(const std::vector<Molecule*> molecules);
    GraphIndexData toIndexData(const std::vector<Assembly*> assemblies);
    graph::Database createGraphData(const GraphObjects& objects);
    assembly::Graph createAssemblyGraph(const assembly::Indices& indices, const graph::Database& atomGraphData);
    assembly::Graph createCompleteAssemblyGraph(const GraphIndexData& data);
    assembly::Graph createVisibleAssemblyGraph(const GraphIndexData& data);
    assembly::Bounds toAssemblyBounds(const assembly::Graph& graph, const std::vector<Sphere>& atomBounds);
} // namespace cds

#endif
