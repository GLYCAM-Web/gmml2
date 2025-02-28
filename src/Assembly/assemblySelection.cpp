#include "includes/Assembly/assemblySelection.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace assembly
{
    Selection selectByAtoms(const Graph& graph, const std::vector<bool>& atoms)
    {
        std::vector<bool> residues = codeUtils::groupContainsSelected(graph.residueCount, graph.atomResidue, atoms);
        std::vector<bool> molecules =
            codeUtils::groupContainsSelected(graph.moleculeCount, graph.residueMolecule, residues);
        return {atoms, residues, molecules};
    }

    Selection selectByMolecules(const Graph& graph, const std::vector<bool>& molecules)
    {
        std::vector<bool> residues = codeUtils::indicesToValues(molecules, graph.residueMolecule);
        std::vector<bool> atoms    = codeUtils::indicesToValues(residues, graph.atomResidue);
        return {atoms, residues, molecules};
    }

    Selection selectByAtomsAndMolecules(const Graph& graph, const std::vector<bool>& atoms,
                                        const std::vector<bool>& molecules)
    {
        Selection byMolecules = selectByMolecules(graph, molecules);
        return selectByAtoms(graph, codeUtils::vectorAnd(atoms, byMolecules.atoms));
    }

    std::vector<size_t> selectedMoleculeResidues(const Graph& graph, const Selection& selection, size_t molecule)
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
} // namespace assembly
