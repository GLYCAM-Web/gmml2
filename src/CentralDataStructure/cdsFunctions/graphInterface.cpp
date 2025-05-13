#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
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
    std::vector<size_t> moleculeAssembly {0};

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

    return {
        {atoms.size(), residues.size(), 1, 1, atomResidue, residueMolecule, moleculeAssembly},
        {atoms, residues, {}, {}}
    };
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
    std::vector<size_t> moleculeAssembly;

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
        moleculeAssembly.push_back(0);
        moleculeIndex++;
    }

    return {
        {atoms.size(), residues.size(), molecules.size(), 1, atomResidue, residueMolecule, moleculeAssembly},
        {atoms, residues, molecules, {}}
    };
}

cds::GraphIndexData cds::toIndexData(const std::vector<Assembly*> assemblies)
{
    size_t assemblyIndex = 0;
    size_t moleculeIndex = 0;
    size_t residueIndex  = 0;
    size_t atomIndex     = 0;

    std::vector<Atom*> atoms;
    std::vector<Residue*> residues;
    std::vector<Molecule*> molecules;
    std::vector<size_t> atomResidue;
    std::vector<size_t> residueMolecule;
    std::vector<size_t> moleculeAssembly;

    for (auto& assembly : assemblies)
    {
        for (auto& molecule : assembly->getMolecules())
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
            moleculeAssembly.push_back(assemblyIndex);
            molecules.push_back(molecule);
            moleculeIndex++;
        }
        assemblyIndex++;
    }

    return {
        {atoms.size(), residues.size(), molecules.size(), assemblies.size(), atomResidue, residueMolecule,
         moleculeAssembly},
        {atoms, residues, molecules, assemblies}
    };
}

graph::Database cds::createGraphData(const GraphObjects& objects)
{
    std::vector<uint> initialIndices;
    // save indices
    for (auto& atom : objects.atoms)
    {
        initialIndices.push_back(atom->getIndex());
    }
    graph::Database graph;
    auto& atoms = objects.atoms;
    for (size_t n = 0; n < atoms.size(); n++)
    {
        atoms[n]->setIndex(n);
        addNode(graph);
    }
    for (size_t n = 0; n < atoms.size(); n++)
    {
        for (auto& neighbor : atoms[n]->getChildren())
        {
            if (codeUtils::contains(objects.atoms, neighbor))
            {
                size_t index = neighbor->getIndex();
                addEdge(graph, {n, index});
            }
        }
    }
    // restore indices
    for (size_t n = 0; n < objects.atoms.size(); n++)
    {
        objects.atoms[n]->setIndex(initialIndices[n]);
    }
    return graph;
}

assembly::Graph cds::createAssemblyGraph(const GraphIndexData& data, const std::vector<bool>& includedAtoms)
{
    const assembly::Indices& indices = data.indices;
    const GraphObjects& objects      = data.objects;
    graph::Database atomGraphData    = cds::createGraphData(objects);
    atomGraphData.nodeAlive          = includedAtoms;
    graph::Graph atomGraph           = graph::identity(atomGraphData);
    graph::Graph residueGraph        = graph::quotient(atomGraphData, indices.atomResidue);
    graph::Database residueData      = graph::asData(residueGraph);
    graph::Graph moleculeGraph =
        graph::quotient(residueData, codeUtils::indicesToValues(indices.residueMolecule, residueData.nodes));
    graph::Database moleculeData = graph::asData(moleculeGraph);
    graph::Graph assemblyGraph =
        graph::quotient(moleculeData, codeUtils::indicesToValues(indices.moleculeAssembly, moleculeData.nodes));
    return assembly::Graph {
        {objects.atoms.size(), objects.residues.size(), objects.molecules.size(), objects.assemblies.size(),
         indices.atomResidue, indices.residueMolecule, indices.moleculeAssembly},
        atomGraph,
        residueGraph,
        moleculeGraph,
        assemblyGraph
    };
}

assembly::Graph cds::createCompleteAssemblyGraph(const GraphIndexData& data)
{
    return createAssemblyGraph(data, std::vector<bool>(data.objects.atoms.size(), true));
}

assembly::Graph cds::createVisibleAssemblyGraph(const GraphIndexData& data)
{
    return createAssemblyGraph(data, cds::atomVisibility(data.objects.atoms));
}

assembly::Bounds cds::toAssemblyBounds(const codeUtils::SparseVector<double>& elementRadii, const GraphIndexData& data,
                                       const assembly::Graph& graph)
{
    std::vector<Sphere> atomBounds = atomCoordinatesWithRadii(elementRadii, data.objects.atoms);
    size_t residueCount            = data.objects.residues.size();
    std::vector<Sphere> residueBounds;
    residueBounds.reserve(residueCount);
    for (size_t n = 0; n < residueCount; n++)
    {
        residueBounds.push_back(boundingSphere(codeUtils::indicesToValues(atomBounds, residueAtoms(graph, n))));
    }
    size_t moleculeCount = data.objects.molecules.size();
    std::vector<Sphere> moleculeBounds;
    moleculeBounds.reserve(moleculeCount);
    for (size_t n = 0; n < moleculeCount; n++)
    {
        moleculeBounds.push_back(boundingSphere(codeUtils::indicesToValues(residueBounds, moleculeResidues(graph, n))));
    }
    return {atomBounds, residueBounds, moleculeBounds};
}
