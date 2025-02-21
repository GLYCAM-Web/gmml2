#ifndef INCLUDES_ASSEMBLY_ASSEMBLYBOUNDS_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYBOUNDS_HPP

#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <vector>

namespace assembly
{
    struct Bounds
    {
        std::vector<cds::Sphere> atoms;
        std::vector<cds::Sphere> residues;
        std::vector<cds::Sphere> molecules;
    };

    void updateResidueBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateResidueMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateMoleculeBounds(const Graph& graph, Bounds& bounds, size_t index);
    void updateBoundsContainingAtoms(const Graph& graph, Bounds& bounds, const std::vector<size_t>& atoms);
} // namespace assembly

#endif