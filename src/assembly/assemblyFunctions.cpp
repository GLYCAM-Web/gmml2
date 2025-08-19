#include "include/assembly/assemblyFunctions.hpp"

#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"

#include <vector>

namespace gmml
{
    namespace assembly
    {
        size_t addAtom(Assembly& assembly, size_t residueId)
        {
            assembly.indices.atomCount++;
            assembly.indices.atomResidue.push_back(residueId);
            assembly.indices.atomAlive.push_back(true);
            return graph::addNode(assembly.atomGraph);
        }

        size_t addResidue(Assembly& assembly, size_t moleculeId)
        {
            size_t id = assembly.indices.residueCount;
            assembly.indices.residueCount++;
            assembly.indices.residueMolecule.push_back(moleculeId);
            return id;
        }

        size_t addMolecule(Assembly& assembly, size_t assemblyId)
        {
            size_t id = assembly.indices.moleculeCount;
            assembly.indices.moleculeCount++;
            assembly.indices.moleculeAssembly.push_back(assemblyId);
            return id;
        }

        size_t addAssembly(Assembly& assembly)
        {
            size_t id = assembly.indices.assemblyCount;
            assembly.indices.assemblyCount++;
            return id;
        }
    } // namespace assembly
} // namespace gmml
