#include "include/assembly/assemblyIndices.hpp"

#include "include/util/containers.hpp"

#include <vector>

namespace gmml
{
    namespace assembly
    {
        const std::vector<size_t>& atomResidues(const Indices& indices) { return indices.atomResidue; }

        const std::vector<size_t>& residueMolecules(const Indices& indices) { return indices.residueMolecule; }

        const std::vector<size_t>& moleculeAssemblies(const Indices& indices) { return indices.moleculeAssembly; }

        std::vector<size_t> indicesOfLivingAtoms(const Indices& indices, const std::vector<bool>& atoms)
        {
            return util::boolsToIndices(util::vectorAnd(indices.atomAlive, atoms));
        }

        std::vector<size_t> atomMolecules(const Indices& indices)
        {
            return util::indicesToValues(residueMolecules(indices), atomResidues(indices));
        }

        std::vector<size_t> atomAssemblies(const Indices& indices)
        {
            return util::indicesToValues(moleculeAssemblies(indices), atomMolecules(indices));
        }

        std::vector<size_t> residueAssemblies(const Indices& indices)
        {
            return util::indicesToValues(moleculeAssemblies(indices), residueMolecules(indices));
        }

        std::vector<bool> isResidueAtom(const Indices& indices, size_t residueId)
        {
            return util::vectorEquals(atomResidues(indices), residueId);
        }

        std::vector<bool> isMoleculeAtom(const Indices& indices, size_t moleculeId)
        {
            return util::vectorEquals(atomMolecules(indices), moleculeId);
        }

        std::vector<bool> isMoleculeResidue(const Indices& indices, size_t moleculeId)
        {
            return util::vectorEquals(residueMolecules(indices), moleculeId);
        }

        std::vector<bool> isAssemblyAtom(const Indices& indices, size_t assemblyId)
        {
            return util::vectorEquals(atomAssemblies(indices), assemblyId);
        }

        std::vector<bool> isAssemblyResidue(const Indices& indices, size_t assemblyId)
        {
            return util::vectorEquals(residueAssemblies(indices), assemblyId);
        }

        std::vector<bool> isAssemblyMolecule(const Indices& indices, size_t assemblyId)
        {
            return util::vectorEquals(moleculeAssemblies(indices), assemblyId);
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
            return util::boolsToIndices(isMoleculeResidue(indices, moleculeId));
        }

        std::vector<size_t> assemblyAtoms(const Indices& indices, size_t assemblyId)
        {
            return indicesOfLivingAtoms(indices, isAssemblyAtom(indices, assemblyId));
        }

        std::vector<size_t> assemblyResidues(const Indices& indices, size_t assemblyId)
        {
            return util::boolsToIndices(isAssemblyResidue(indices, assemblyId));
        }

        std::vector<size_t> assemblyMolecules(const Indices& indices, size_t assemblyId)
        {
            return util::boolsToIndices(isAssemblyMolecule(indices, assemblyId));
        }
    } // namespace assembly
} // namespace gmml
