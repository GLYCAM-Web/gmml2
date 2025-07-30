#ifndef INCLUDE_ASSEMBLY_ASSEMBLYBOUNDS_HPP
#define INCLUDE_ASSEMBLY_ASSEMBLYBOUNDS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <vector>

namespace gmml
{
    namespace assembly
    {
        void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& selectedAtoms);
        Bounds toAssemblyBounds(const Graph& graph, const std::vector<Sphere>& atomBounds);
    } // namespace assembly
} // namespace gmml

#endif