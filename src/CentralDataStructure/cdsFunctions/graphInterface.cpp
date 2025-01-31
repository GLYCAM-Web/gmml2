#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

cds::GraphIndexData cds::toIndexData(const std::vector<Residue*> inputResidues)
{
    size_t residueIndex = 0;
    size_t atomIndex    = 0;

    std::vector<Atom*> atoms;
    std::vector<Residue*> residues;
    std::vector<size_t> atomResidue;
    std::vector<size_t> residueMolecule;

    for (auto& residue : inputResidues)
    {
        for (auto& atom : residue->getAtoms())
        {
            atomResidue.push_back(residueIndex);
            atoms.push_back(atom);
            atomIndex++;
        }
        residueMolecule.push_back(0);
        residues.push_back(residue);
        residueIndex++;
    }

    return {atoms, residues, {}, atomResidue, residueMolecule};
}

cds::GraphIndexData cds::toIndexData(const std::vector<Molecule*> molecules)
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

graph::Database cds::createGraphData(const GraphIndexData& indices)
{
    std::vector<uint> initialIndices;
    // save indices
    for (auto& atom : indices.atoms)
    {
        initialIndices.push_back(atom->getIndex());
    }
    graph::Database graph;
    auto& atoms = indices.atoms;
    for (size_t n = 0; n < atoms.size(); n++)
    {
        atoms[n]->setIndex(n);
        addNode(graph);
    }
    for (size_t n = 0; n < atoms.size(); n++)
    {
        for (auto& neighbor : atoms[n]->getChildren())
        {
            if (codeUtils::contains(indices.atoms, neighbor))
            {
                size_t index = neighbor->getIndex();
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

assembly::Graph cds::createAssemblyGraph(const GraphIndexData& indices)
{
    graph::Database atomGraphData = cds::createGraphData(indices);
    graph::Graph atomGraph        = graph::identity(atomGraphData);
    graph::Graph residueGraph     = graph::quotient(atomGraphData, indices.atomResidue);
    graph::Graph moleculeGraph    = graph::quotient(graph::asData(residueGraph), indices.residueMolecule);
    return assembly::Graph {indices.atoms.size(),
                            indices.residues.size(),
                            indices.molecules.size(),
                            indices.atomResidue,
                            indices.residueMolecule,
                            atomGraph,
                            residueGraph,
                            moleculeGraph};
}
