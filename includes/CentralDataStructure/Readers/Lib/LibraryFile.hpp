#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYFILE_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <istream>
#include <string>
#include <array>
#include <vector>

namespace lib
{
    struct AtomData
    {
        std::vector<size_t> residueIndex;
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<double> charges;
        std::vector<int> numbers;
        std::vector<cds::Coordinate> coordinates;
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<bool> hasCoordinates;
    };

    struct LibraryData
    {
        AtomData atoms;
        ResidueData residues;
        std::vector<std::array<size_t, 2>> bonds;
    };

    void parseMolecule(cds::Molecule* molecule, const std::string& filename);

} // namespace lib
#endif
