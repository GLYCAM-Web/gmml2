#ifndef INCLUDES_ASSEMBLY_ASSEMBLYINDICES_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYINDICES_HPP

#include <cstddef>
#include <vector>

namespace assembly
{
    struct Indices
    {
        size_t atomCount     = 0;
        size_t residueCount  = 0;
        size_t moleculeCount = 0;
        size_t assemblyCount = 0;
        std::vector<bool> atomAlive;
        std::vector<size_t> atomResidue;
        std::vector<size_t> residueMolecule;
        std::vector<size_t> moleculeAssembly;
    };
} // namespace assembly

#endif
