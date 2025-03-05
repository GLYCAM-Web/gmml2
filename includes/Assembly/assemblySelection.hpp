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

    Selection selectAll(const Graph& graph);
    Selection selectByAtoms(const Graph& graph, const std::vector<bool>& atoms);
    Selection selectByResidues(const Graph& graph, const std::vector<bool>& residues);
    Selection selectByMolecules(const Graph& graph, const std::vector<bool>& molecules);
    Selection selectByAtomsAndMolecules(const Graph& graph, const std::vector<bool>& atoms,
                                        const std::vector<bool>& molecules);

    Selection intersection(const Graph& graph, const Selection& first, const Selection& second);

    std::vector<size_t> selectedMolecules(const Selection& selection);
    std::vector<size_t> moleculeSelectedResidues(const Graph& graph, const Selection& selection, size_t molecule);
    std::vector<std::vector<size_t>> moleculeSelectedResidues(const Graph& graph, const Selection& selection,
                                                              const std::vector<size_t>& molecules);
    std::vector<size_t> moleculeSelectedAtoms(const Graph& graph, const Selection& selection, size_t molecule);
} // namespace assembly

#endif