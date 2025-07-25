#ifndef INCLUDE_READERS_LIB_LIBRARYFILE_HPP
#define INCLUDE_READERS_LIB_LIBRARYFILE_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"

#include <array>
#include <string>
#include <vector>

namespace gmml
{
    namespace lib
    {
        struct AtomData
        {
            std::vector<std::string> names;
            std::vector<std::string> types;
            std::vector<Element> elements;
            std::vector<double> charges;
            std::vector<uint> numbers;
            std::vector<Coordinate> coordinates;
        };

        struct ResidueData
        {
            bool hasCoordinates;
            AtomData atoms;
            std::vector<std::array<size_t, 2>> bonds;
        };

        struct LibraryData
        {
            std::vector<std::string> residueNames;
            std::vector<ResidueData> residues;
        };

        LibraryData loadLibraryData(const std::string& filename);

    } // namespace lib
} // namespace gmml

#endif
