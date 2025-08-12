#ifndef INCLUDE_CARBOHYDRATE_CARBOHYDRATETYPES_HPP
#define INCLUDE_CARBOHYDRATE_CARBOHYDRATETYPES_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/elements.hpp"
#include "include/metadata/residueTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace carbohydrate
    {
        struct AtomData
        {
            std::vector<std::string> names;
            std::vector<std::string> types;
            std::vector<uint> numbers;
            std::vector<uint> atomicNumbers;
            std::vector<Element> elements;
            std::vector<Coordinate> coordinates;
            std::vector<double> charges;
            std::vector<bool> visible;
        };

        struct ResidueData
        {
            std::vector<std::string> names;
            std::vector<ResidueType> types;
            std::vector<std::string> ids;
            std::vector<uint> numbers;
        };

        struct EdgeData
        {
            std::vector<std::string> names;
        };

        struct CarbohydrateData
        {
            AtomData atoms;
            ResidueData residues;
            EdgeData edges;
            assembly::Indices indices;
            graph::Database atomGraph;
        };
    } // namespace carbohydrate
} // namespace gmml

#endif
