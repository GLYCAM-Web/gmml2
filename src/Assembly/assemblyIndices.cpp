#include "includes/Assembly/assemblyIndices.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace assembly
{
    const std::vector<size_t>& atomResidues(const Indices& indices)
    {
        return indices.atomResidue;
    }

    const std::vector<size_t>& residueMolecules(const Indices& indices)
    {
        return indices.residueMolecule;
    }

    const std::vector<size_t>& moleculeAssemblies(const Indices& indices)
    {
        return indices.moleculeAssembly;
    }

    std::vector<size_t> indicesOfLivingAtoms(const Indices& indices, const std::vector<bool>& atoms)
    {
        return codeUtils::boolsToIndices(codeUtils::vectorAnd(indices.atomAlive, atoms));
    }

    std::vector<size_t> atomMolecules(const Indices& indices)
    {
        return codeUtils::indicesToValues(residueMolecules(indices), atomResidues(indices));
    }

    std::vector<size_t> atomAssemblies(const Indices& indices)
    {
        return codeUtils::indicesToValues(moleculeAssemblies(indices), atomMolecules(indices));
    }

    std::vector<size_t> residueAssemblies(const Indices& indices)
    {
        return codeUtils::indicesToValues(moleculeAssemblies(indices), residueMolecules(indices));
    }

    std::vector<bool> isResidueAtom(const Indices& indices, size_t residueId)
    {
        return codeUtils::vectorEquals(atomResidues(indices), residueId);
    }

    std::vector<bool> isMoleculeAtom(const Indices& indices, size_t moleculeId)
    {
        return codeUtils::vectorEquals(atomMolecules(indices), moleculeId);
    }

    std::vector<bool> isMoleculeResidue(const Indices& indices, size_t moleculeId)
    {
        return codeUtils::vectorEquals(residueMolecules(indices), moleculeId);
    }

    std::vector<bool> isAssemblyAtom(const Indices& indices, size_t assemblyId)
    {
        return codeUtils::vectorEquals(atomAssemblies(indices), assemblyId);
    }

    std::vector<bool> isAssemblyResidue(const Indices& indices, size_t assemblyId)
    {
        return codeUtils::vectorEquals(residueAssemblies(indices), assemblyId);
    }

    std::vector<bool> isAssemblyMolecule(const Indices& indices, size_t assemblyId)
    {
        return codeUtils::vectorEquals(moleculeAssemblies(indices), assemblyId);
    }

    std::vector<size_t> residueAtoms(const Indices& indices, size_t residueId)
    {
        return indicesOfLivingAtoms(indices, isResidueAtom(indices, residueId));
    }

    std::vector<size_t> moleculeAtoms(const Indices& indices, size_t moleculeId)
    {
        return indicesOfLivingAtoms(indices, isMoleculeAtom(indices, moleculeId));
    }

    std::vector<size_t> moleculeResidues(const Indices& indices, size_t moleculeId)
    {
        return codeUtils::boolsToIndices(isMoleculeResidue(indices, moleculeId));
    }

    std::vector<size_t> assemblyAtoms(const Indices& indices, size_t assemblyId)
    {
        return indicesOfLivingAtoms(indices, isAssemblyAtom(indices, assemblyId));
    }

    std::vector<size_t> assemblyResidues(const Indices& indices, size_t assemblyId)
    {
        return codeUtils::boolsToIndices(isAssemblyResidue(indices, assemblyId));
    }

    std::vector<size_t> assemblyMolecules(const Indices& indices, size_t assemblyId)
    {
        return codeUtils::boolsToIndices(isAssemblyMolecule(indices, assemblyId));
    }

} // namespace assembly