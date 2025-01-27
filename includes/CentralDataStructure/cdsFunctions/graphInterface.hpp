#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_GRAPHINTERFACE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <vector>

namespace cds
{
    struct GraphIndexData
    {
        GraphIndexData(const std::vector<Atom*>& atoms_, const std::vector<Residue*>& residues_,
                       const std::vector<Molecule*>& molecules_, const std::vector<size_t>& atomResidue_,
                       const std::vector<size_t>& residueMolecule_)
            : atoms(atoms_), residues(residues_), molecules(molecules_), atomResidue(atomResidue_),
              residueMolecule(residueMolecule_)
        {}

        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
    };

    GraphIndexData toIndexData(std::vector<Molecule*> molecules);
    graph::Database createGraphData(GraphIndexData& indices);
    assembly::Graph createAssemblyGraph(GraphIndexData& indices);
} // namespace cds

#endif
