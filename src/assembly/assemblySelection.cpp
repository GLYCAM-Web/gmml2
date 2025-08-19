#include "include/assembly/assemblySelection.hpp"

#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/util/containers.hpp"

#include <vector>

namespace gmml
{
    namespace assembly
    {
        Selection selectAll(const Graph& graph)
        {
            std::vector<bool> atoms(atomCount(graph.source), true);
            std::vector<bool> residues(residueCount(graph.source), true);
            std::vector<bool> molecules(moleculeCount(graph.source), true);
            return {atoms, residues, molecules};
        }

        Selection selectByAtoms(const Graph& graph, const std::vector<bool>& atoms)
        {
            Selection all = selectAll(graph);
            Selection some = all;
            some.atoms = atoms;
            return intersection(graph, some, all);
        }

        Selection selectByResidues(const Graph& graph, const std::vector<bool>& residues)
        {
            Selection all = selectAll(graph);
            Selection some = all;
            some.residues = residues;
            return intersection(graph, some, all);
        }

        Selection selectByMolecules(const Graph& graph, const std::vector<bool>& molecules)
        {
            Selection all = selectAll(graph);
            Selection some = all;
            some.molecules = molecules;
            return intersection(graph, some, all);
        }

        Selection selectByAtomsAndMolecules(
            const Graph& graph, const std::vector<bool>& atoms, const std::vector<bool>& molecules)
        {
            Selection all = selectAll(graph);
            Selection some = all;
            some.atoms = atoms;
            some.molecules = molecules;
            return intersection(graph, some, all);
        }

        Selection intersection(const Graph& graph, const Selection& first, const Selection& second)
        {
            auto setValues =
                [](std::vector<bool>& result, const std::vector<size_t>& group, const std::vector<bool>& value)
            {
                std::fill(result.begin(), result.end(), false);
                for (size_t n = 0; n < value.size(); n++)
                {
                    size_t index = group[n];
                    result[index] = result[index] || value[n];
                }
            };
            size_t atomCount = assembly::atomCount(graph.source);
            size_t residueCount = assembly::residueCount(graph.source);
            size_t moleculeCount = assembly::moleculeCount(graph.source);
            const Indices& indices = graph.source.indices;
            std::vector<bool> atoms(atomCount, false);
            std::vector<bool> residues(residueCount, false);
            std::vector<bool> molecules(moleculeCount, false);

            for (size_t n = 0; n < moleculeCount; n++)
            {
                molecules[n] = first.molecules[n] && second.molecules[n];
            }
            for (size_t n = 0; n < residueCount; n++)
            {
                residues[n] = molecules[residueMolecule(indices, n)] && first.residues[n] && second.residues[n];
            }
            for (size_t n = 0; n < atomCount; n++)
            {
                atoms[n] = residues[atomResidue(indices, n)] && first.atoms[n] && second.atoms[n];
            }
            setValues(residues, indices.atomResidue, atoms);
            setValues(molecules, indices.residueMolecule, residues);
            return {atoms, residues, molecules};
        }

        std::vector<size_t> selectedMolecules(const Selection& selection)
        {
            return util::boolsToIndices(selection.molecules);
        }

        std::vector<size_t> moleculeSelectedResidues(const Graph& graph, const Selection& selection, size_t molecule)
        {
            const std::vector<size_t>& residues = moleculeResidues(graph, molecule);
            std::vector<size_t> selectedResidues;
            selectedResidues.reserve(residues.size());
            for (size_t n : residues)
            {
                if (selection.residues[n])
                {
                    selectedResidues.push_back(n);
                }
            }
            return selectedResidues;
        }

        std::vector<std::vector<size_t>> moleculeSelectedResidues(
            const Graph& graph, const Selection& selection, const std::vector<size_t>& molecules)
        {
            std::vector<std::vector<size_t>> result;
            result.reserve(molecules.size());
            for (size_t molecule : molecules)
            {
                result.push_back(moleculeSelectedResidues(graph, selection, molecule));
            }
            return result;
        }

        std::vector<size_t> moleculeSelectedAtoms(const Graph& graph, const Selection& selection, size_t molecule)
        {
            const std::vector<size_t>& atoms = moleculeAtoms(graph, molecule);
            std::vector<size_t> selectedAtoms;
            selectedAtoms.reserve(atoms.size());
            for (size_t n : atoms)
            {
                if (selection.atoms[n])
                {
                    selectedAtoms.push_back(n);
                }
            }
            return selectedAtoms;
        }
    } // namespace assembly
} // namespace gmml
