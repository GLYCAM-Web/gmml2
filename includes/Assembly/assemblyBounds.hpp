#ifndef INCLUDES_ASSEMBLY_ASSEMBLYBOUNDS_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYBOUNDS_HPP

#include "includes/Assembly/assemblyTypes.hpp"

#include <vector>

namespace assembly
{
    void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& selectedAtoms);
} // namespace assembly

#endif