#ifndef INCLUDES_ASSEMBLY_ASSEMBLYTYPES_HPP
#define INCLUDES_ASSEMBLY_ASSEMBLYTYPES_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"

#include <array>
#include <vector>

namespace gmml
{
    namespace assembly
    {
        struct Indices
        {
            size_t atomCount = 0;
            size_t residueCount = 0;
            size_t moleculeCount = 0;
            size_t assemblyCount = 0;
            std::vector<bool> atomAlive;
            std::vector<size_t> atomResidue;
            std::vector<size_t> residueMolecule;
            std::vector<size_t> moleculeAssembly;
        };

        struct Bounds
        {
            std::vector<Sphere> atoms;
            std::vector<Sphere> residues;
            std::vector<Sphere> molecules;
        };

        struct Selection
        {
            std::vector<bool> atoms;
            std::vector<bool> residues;
            std::vector<bool> molecules;
        };

        struct Graph
        {
            Indices indices;
            graph::Graph atoms;
            graph::Graph residues;
            graph::Graph molecules;
            graph::Graph assemblies;
        };
    } // namespace assembly
} // namespace gmml

#endif
