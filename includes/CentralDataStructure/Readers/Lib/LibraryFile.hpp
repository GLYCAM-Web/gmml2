#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/MolecularMetadata/elements.hpp"

#include <istream>
#include <string>
#include <array>
#include <vector>

namespace lib
{
    struct AtomData
    {
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<MolecularMetadata::Element> elements;
        std::vector<double> charges;
        std::vector<uint> numbers;
        std::vector<cds::Coordinate> coordinates;
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
#endif
