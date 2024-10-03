#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/Graph/types.hpp"
#include "includes/Graph/manipulation.hpp"

#include <vector>

cds::GraphIndexData cds::toIndexData(std::vector<Molecule*> molecules)
{
    size_t moleculeIndex = 0;
    size_t residueIndex  = 0;
    size_t atomIndex     = 0;

    std::vector<Atom*> atoms;
    std::vector<Residue*> residues;
    std::vector<size_t> atomResidue;
    std::vector<size_t> residueMolecule;

    for (auto& molecule : molecules)
    {
        for (auto& residue : molecule->getResidues())
        {
            for (auto& atom : residue->getAtoms())
            {
                atomResidue.push_back(residueIndex);
                atoms.push_back(atom);
                atomIndex++;
            }
            residueMolecule.push_back(moleculeIndex);
            residues.push_back(residue);
            residueIndex++;
        }
        moleculeIndex++;
    }

    return {atoms, residues, molecules, atomResidue, residueMolecule};
}

graph::Database cds::createGraphData(GraphIndexData& indices)
{
    std::vector<uint> initialIndices;
    // save indices
    for (auto& atom : indices.atoms)
    {
        initialIndices.push_back(atom->getIndex());
    }
    graph::Database graph({}, {}, {}, {}, {});
    auto& atoms = indices.atoms;
    for (size_t n = 0; n < atoms.size(); n++)
    {
        atoms[n]->setIndex(n);
        addNode(graph);
    }
    for (size_t n = 0; n < atoms.size(); n++)
    {
        for (auto& neighbor : atoms[n]->getNeighbors())
        {
            size_t index = neighbor->getIndex();
            if (n < index)
            {
                addEdge(graph, {n, index});
            }
        }
    }
    // restore indices
    for (size_t n = 0; n < indices.atoms.size(); n++)
    {
        indices.atoms[n]->setIndex(initialIndices[n]);
    }
    return graph;
}
