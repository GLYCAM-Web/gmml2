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
    namespace
    {
        void sortByIndex(std::vector<Atom*>& atoms, std::vector<Residue*>& residues, std::vector<Molecule*>& molecules)
        {
            std::sort(
                atoms.begin(), atoms.end(), [](Atom* a, Atom* b) { return a->getDataIndex() < b->getDataIndex(); });
            std::sort(
                residues.begin(),
                residues.end(),
                [](Residue* a, Residue* b) { return a->getDataIndex() < b->getDataIndex(); });
            std::sort(
                molecules.begin(),
                molecules.end(),
                [](Molecule* a, Molecule* b) { return a->getDataIndex() < b->getDataIndex(); });
        }
    } // namespace

    AssemblyIndexOffset reorderDataIndices(std::vector<Molecule*>& molecules, AssemblyIndexOffset offset)
    {
        for (auto molecule : molecules)
        {
            molecule->setDataIndex(offset.molecule);
            offset.molecule++;
            for (auto residue : molecule->getResidues())
            {
                residue->setDataIndex(offset.residue);
                offset.residue++;
                for (auto atom : residue->getAtoms())
                {
                    atom->setDataIndex(offset.atom);
                    offset.atom++;
                }
            }
        }
        return offset;
    }

    GraphIndexData toIndexData(const std::vector<Residue*> inputResidues)
    {
        size_t residueIndex = 0;
        size_t atomIndex = 0;

        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly {0};

        residues.reserve(inputResidues.size());
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

        sortByIndex(atoms, residues, molecules);

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

    GraphIndexData toIndexData(const std::vector<Molecule*> inputMolecules)
    {
        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly;

        molecules.reserve(inputMolecules.size());
        for (auto& molecule : inputMolecules)
        {
            size_t moleculeIndex = molecule->getDataIndex();
            for (auto& residue : molecule->getResidues())
            {
                size_t residueIndex = residue->getDataIndex();
                for (auto& atom : residue->getAtoms())
                {
                    atomResidue.push_back(residueIndex);
                    atoms.push_back(atom);
                }
                residueMolecule.push_back(moleculeIndex);
                residues.push_back(residue);
            }
            moleculeAssembly.push_back(0);
            molecules.push_back(molecule);
        }

        sortByIndex(atoms, residues, molecules);

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
                size_t moleculeIndex = molecule->getDataIndex();
                for (auto& residue : molecule->getResidues())
                {
                    size_t residueIndex = residue->getDataIndex();
                    for (auto& atom : residue->getAtoms())
                    {
                        atomResidue.push_back(residueIndex);
                        atoms.push_back(atom);
                    }
                    residueMolecule.push_back(moleculeIndex);
                    residues.push_back(residue);
                }
                moleculeAssembly.push_back(assemblyIndex);
                molecules.push_back(molecule);
            }
            assemblyIndex++;
        }

        sortByIndex(atoms, residues, molecules);

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
            initialIndices.push_back(atom->getDataIndex());
        }
        graph::Database graph;
        for (size_t n = 0; n < atoms.size(); n++)
        {
            atoms[n]->setIndex(n);
            addNode(graph);
            graph.nodes.alive[n] = !atoms[n]->isSoftDeleted();
        }
        for (size_t n = 0; n < atoms.size(); n++)
        {
            for (auto& neighbor : atoms[n]->getChildren())
            {
                if (util::contains(atoms, neighbor))
                {
                    size_t index = neighbor->getDataIndex();
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

    assembly::Graph createCompleteAssemblyGraph(const GraphIndexData& data)
    {
        return createAssemblyGraph(data.indices, createGraphData(data.objects));
    }
} // namespace gmml
