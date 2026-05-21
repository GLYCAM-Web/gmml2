#ifndef INCLUDE_CENTRALDATASTRUCTURE_GRAPHINTERFACE_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_GRAPHINTERFACE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/util/containerTypes.hpp"

#include <vector>

namespace gmml
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

    struct AssemblyIndexOffset
    {
        size_t molecule;
        size_t residue;
        size_t atom;
    };

    AssemblyIndexOffset reorderDataIndices(std::vector<Molecule*>& molecules, AssemblyIndexOffset offset);
    GraphIndexData toIndexData(const std::vector<Residue*> inputResidues);
    GraphIndexData toIndexData(const std::vector<Molecule*> molecules);
    GraphIndexData toIndexData(const std::vector<Assembly*> assemblies);
    graph::Database createGraphData(const GraphObjects& objects);
    assembly::Graph createCompleteAssemblyGraph(const GraphIndexData& data);
} // namespace gmml

#endif
