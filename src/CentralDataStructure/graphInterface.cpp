#include "include/CentralDataStructure/graphInterface.hpp"

#include "include/CentralDataStructure/assembly.hpp"
#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"

#include <vector>

namespace gmml
{
    GraphIndexData toIndexData(const std::vector<Residue*> inputResidues)
    {
        size_t residueIndex = 0;
        size_t atomIndex = 0;

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
            {atoms.size(),
             residues.size(),
             1, 1,
             std::vector<bool>(true, atoms.size()),
             atomResidue, residueMolecule,
             moleculeAssembly},
            {atoms, residues, {}, {}}
        };
    }

    GraphIndexData toIndexData(const std::vector<Molecule*> molecules)
    {
        size_t moleculeIndex = 0;
        size_t residueIndex = 0;
        size_t atomIndex = 0;

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
            {atoms.size(),
             residues.size(),
             molecules.size(),
             1, std::vector<bool>(atoms.size(), true),
             atomResidue, residueMolecule,
             moleculeAssembly},
            {atoms, residues, molecules, {}}
        };
    }

    GraphIndexData toIndexData(const std::vector<Assembly*> assemblies)
    {
        size_t assemblyIndex = 0;
        size_t moleculeIndex = 0;
        size_t residueIndex = 0;
        size_t atomIndex = 0;

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
            {atoms.size(),
             residues.size(),
             molecules.size(),
             assemblies.size(),
             std::vector<bool>(atoms.size(), true),
             atomResidue, residueMolecule,
             moleculeAssembly},
            {atoms, residues, molecules, assemblies}
        };
    }

    graph::Database createGraphData(const GraphObjects& objects)
    {
        std::vector<uint> initialIndices;
        const std::vector<Atom*>& atoms = objects.atoms;
        // save indices
        for (auto& atom : atoms)
        {
            initialIndices.push_back(atom->getIndex());
        }
        graph::Database graph;
        for (size_t n = 0; n < atoms.size(); n++)
        {
            atoms[n]->setIndex(n);
            addNode(graph);
        }
        for (size_t n = 0; n < atoms.size(); n++)
        {
            for (auto& neighbor : atoms[n]->getChildren())
            {
                if (util::contains(atoms, neighbor))
                {
                    size_t index = neighbor->getIndex();
                    addEdge(graph, {n, index});
                }
            }
        }
        // restore indices
        for (size_t n = 0; n < atoms.size(); n++)
        {
            atoms[n]->setIndex(initialIndices[n]);
        }
        return graph;
    }

    assembly::Graph createAssemblyGraph(const assembly::Indices& indices, const graph::Database& atomGraphData)
    {
        graph::Graph atomGraph = graph::identity(atomGraphData);
        graph::Graph residueGraph = graph::quotient(atomGraphData, indices.atomResidue);
        graph::Database residueData = graph::asData(residueGraph);
        graph::Graph moleculeGraph =
            graph::quotient(residueData, util::indicesToValues(indices.residueMolecule, residueData.nodes));
        graph::Database moleculeData = graph::asData(moleculeGraph);
        graph::Graph assemblyGraph =
            graph::quotient(moleculeData, util::indicesToValues(indices.moleculeAssembly, moleculeData.nodes));
        return assembly::Graph {
            {indices.atomCount,
             indices.residueCount,
             indices.moleculeCount,
             indices.assemblyCount,
             atomGraphData.nodeAlive,
             indices.atomResidue,
             indices.residueMolecule,
             indices.moleculeAssembly},
            atomGraph,
            residueGraph,
            moleculeGraph,
            assemblyGraph
        };
    }

    assembly::Graph createCompleteAssemblyGraph(const GraphIndexData& data)
    {
        return createAssemblyGraph(data.indices, createGraphData(data.objects));
    }

    assembly::Graph createVisibleAssemblyGraph(const GraphIndexData& data)
    {
        graph::Database atomGraphData = createGraphData(data.objects);
        atomGraphData.nodeAlive = atomVisibility(data.objects.atoms);
        return createAssemblyGraph(data.indices, atomGraphData);
    }

    assembly::Bounds toAssemblyBounds(const assembly::Graph& graph, const std::vector<Sphere>& atomBounds)
    {
        size_t residueCount = graph.indices.residueCount;
        std::vector<Sphere> residueBounds;
        residueBounds.reserve(residueCount);
        for (size_t n = 0; n < residueCount; n++)
        {
            residueBounds.push_back(boundingSphere(util::indicesToValues(atomBounds, residueAtoms(graph, n))));
        }
        size_t moleculeCount = graph.indices.moleculeCount;
        std::vector<Sphere> moleculeBounds;
        moleculeBounds.reserve(moleculeCount);
        for (size_t n = 0; n < moleculeCount; n++)
        {
            moleculeBounds.push_back(boundingSphere(util::indicesToValues(residueBounds, moleculeResidues(graph, n))));
        }
        return {atomBounds, residueBounds, moleculeBounds};
    }
} // namespace gmml
