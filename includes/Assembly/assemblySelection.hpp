#ifndef INCLUDES_ASSEMBLY_ASSEMBLYSELECTION_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYSELECTION_HPP

#include "includes/Assembly/assemblyGraph.hpp"

#include <vector>

namespace assembly
{
    struct Selection
    {
        std::vector<bool> atoms;
        std::vector<bool> residues;
        std::vector<bool> molecules;
    };

    Selection selectByAtoms(const Graph& graph, const std::vector<bool>& atoms);
    Selection selectByMolecules(const Graph& graph, const std::vector<bool>& molecules);
    Selection selectByAtomsAndMolecules(const Graph& graph, const std::vector<bool>& atoms,
                                        const std::vector<bool>& molecules);

    std::vector<size_t> selectedMoleculeResidues(const Graph& graph, const Selection& selection, size_t molecule);
} // namespace assembly

#endif