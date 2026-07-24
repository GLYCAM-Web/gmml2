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

        void reindex(std::vector<Atom*>& atoms, std::vector<Residue*>& residues, std::vector<Molecule*>& molecules)
        {
            for (size_t n = 0; n < atoms.size(); n++)
            {
                atoms[n]->setDataIndex(n);
            }
            for (size_t n = 0; n < residues.size(); n++)
            {
                residues[n]->setDataIndex(n);
            }
            for (size_t n = 0; n < molecules.size(); n++)
            {
                molecules[n]->setDataIndex(n);
            }
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
        std::vector<Atom*> atoms;
        std::vector<Residue*> residues = inputResidues;
        std::vector<Molecule*> molecules {};

        for (auto& residue : residues)
        {
            for (auto& atom : residue->getAtoms())
            {
                atoms.push_back(atom);
            }
        }

        sortByIndex(atoms, residues, molecules);
        reindex(atoms, residues, molecules);

        std::vector<size_t> atomResidue(atoms.size(), residues.size());
        std::vector<size_t> residueMolecule(residues.size(), 0);
        std::vector<size_t> moleculeAssembly {0};

        residues.reserve(inputResidues.size());
        for (auto& residue : inputResidues)
        {
            size_t residueIndex = residue->getDataIndex();
            for (auto& atom : residue->getAtoms())
            {
                atomResidue[atom->getDataIndex()] = residueIndex;
            }
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

    GraphIndexData toIndexData(const std::vector<Molecule*> inputMolecules)
    {
        std::vector<Atom*> atoms;
        std::vector<Residue*> residues;
        std::vector<Molecule*> molecules = inputMolecules;

        for (auto& molecule : molecules)
        {
            for (auto& residue : molecule->getResidues())
            {
                for (auto& atom : residue->getAtoms())
                {
                    atoms.push_back(atom);
                }
                residues.push_back(residue);
            }
        }

        sortByIndex(atoms, residues, molecules);
        reindex(atoms, residues, molecules);

        std::vector<size_t> atomResidue(atoms.size(), residues.size());
        std::vector<size_t> residueMolecule(residues.size(), molecules.size());
        std::vector<size_t> moleculeAssembly(molecules.size(), 0);

        for (auto& molecule : molecules)
        {
            size_t moleculeIndex = molecule->getDataIndex();
            for (auto& residue : molecule->getResidues())
            {
                size_t residueIndex = residue->getDataIndex();
                for (auto& atom : residue->getAtoms())
                {
                    atomResidue[atom->getDataIndex()] = residueIndex;
                }
                residueMolecule[residueIndex] = moleculeIndex;
            }
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
        reindex(atoms, residues, molecules);

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
        const std::vector<Atom*>& atoms = objects.atoms;
        graph::Database graph;
        for (size_t n = 0; n < atoms.size(); n++)
        {
            atoms[n]->setDataIndex(n);
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
        return graph;
    }

    assembly::Graph createCompleteAssemblyGraph(const GraphIndexData& data)
    {
        return createAssemblyGraph(data.indices, createGraphData(data.objects));
    }
} // namespace gmml
