#ifndef INCLUDE_CARBOHYDRATE_CARBOHYDRATETYPES_HPP
#define INCLUDE_CARBOHYDRATE_CARBOHYDRATETYPES_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"
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

        struct RotatableBondData
        {
            std::vector<std::array<size_t, 4>> dihedralAtoms;
            std::vector<std::vector<size_t>> movingAtoms;
            std::vector<std::vector<size_t>> metadata;
        };

        struct ResidueLinkageData
        {
            std::vector<RotamerType> rotamerTypes;
            std::vector<size_t> edgeId;
            std::vector<std::vector<size_t>> rotatableBonds;
        };

        struct CarbohydrateData
        {
            AtomData atoms;
            ResidueData residues;
            EdgeData edges;
            RotatableBondData rotatableBonds;
            ResidueLinkageData residueLinkages;
            assembly::Indices indices;
            graph::Database atomGraph;
        };
    } // namespace carbohydrate
} // namespace gmml

#endif
