#ifndef INCLUDE_ASSEMBLY_ASSEMBLYFUNCTIONS_HPP
#define INCLUDE_ASSEMBLY_ASSEMBLYFUNCTIONS_HPP

#include "include/assembly/assemblyTypes.hpp"

namespace gmml
{
    namespace assembly
    {
        size_t addAtom(Assembly& assembly, size_t residueId);
        size_t addResidue(Assembly& assembly, size_t moleculeId);
        size_t addMolecule(Assembly& assembly, size_t assemblyId);
        size_t addAssembly(Assembly& assembly);
    } // namespace assembly
} // namespace gmml

#endif
