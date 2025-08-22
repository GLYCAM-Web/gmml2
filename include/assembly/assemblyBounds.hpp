#ifndef INCLUDE_ASSEMBLY_ASSEMBLYBOUNDS_HPP
#define INCLUDE_ASSEMBLY_ASSEMBLYBOUNDS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/containerTypes.hpp"

#include <vector>

namespace gmml
{
    namespace assembly
    {
        void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
        void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& selectedAtoms);
        std::vector<Sphere> toAtomBounds(
            const util::SparseVector<double>& elementRadii,
            const std::vector<Element>& elements,
            const std::vector<Coordinate>& coordinates);
        Bounds toAssemblyBounds(const Indices& indices, const std::vector<Sphere>& atomBounds);
    } // namespace assembly
} // namespace gmml

#endif