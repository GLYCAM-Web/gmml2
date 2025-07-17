#ifndef INCLUDES_ASSEMBLY_ASSEMBLYINDICES_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYINDICES_HPP

#include "includes/Assembly/assemblyTypes.hpp"

#include <cstddef>
#include <vector>

namespace assembly
{
    const std::vector<size_t>& atomResidues(const Indices& indices);
    const std::vector<size_t>& residueMolecules(const Indices& indices);
    const std::vector<size_t>& moleculeAssemblies(const Indices& indices);
    std::vector<size_t> indicesOfLivingAtoms(const Indices& indices, const std::vector<bool>& atoms);
    std::vector<size_t> atomMolecules(const Indices& indices);
    std::vector<size_t> atomAssemblies(const Indices& indices);
    std::vector<size_t> residueAssemblies(const Indices& indices);
    std::vector<bool> isResidueAtom(const Indices& indices, size_t residueId);
    std::vector<bool> isMoleculeAtom(const Indices& indices, size_t residueId);
    std::vector<bool> isMoleculeResidue(const Indices& indices, size_t residueId);
    std::vector<bool> isAssemblyAtom(const Indices& indices, size_t residueId);
    std::vector<bool> isAssemblyResidue(const Indices& indices, size_t assemblyId);
    std::vector<bool> isAssemblyMolecule(const Indices& indices, size_t assemblyId);
    std::vector<size_t> residueAtoms(const Indices& indices, size_t residueId);
    std::vector<size_t> moleculeAtoms(const Indices& indices, size_t moleculeId);
    std::vector<size_t> moleculeResidues(const Indices& indices, size_t moleculeId);
    std::vector<size_t> assemblyAtoms(const Indices& indices, size_t assemblyId);
    std::vector<size_t> assemblyResidues(const Indices& indices, size_t assemblyId);
    std::vector<size_t> assemblyMolecules(const Indices& indices, size_t assemblyId);
} // namespace assembly

#endif
