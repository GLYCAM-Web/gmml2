#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEDATA_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/residueTypes.hpp"

#include <vector>
#include <string>

namespace cds
{
    struct OffFileAtomData
    {
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<int> atomicNumbers;
        std::vector<double> charges;
        std::vector<Coordinate> coordinates;
        std::vector<size_t> residues;
        std::vector<std::pair<size_t, size_t>> bonds;
    };

    struct OffFileResidueData
    {
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<ResidueType> types;
        std::vector<std::vector<size_t>> atomIndices;
        std::vector<std::vector<size_t>> atomsConnectedToOtherResidues;
    };

    struct OffFileData
    {
        OffFileResidueData residues;
        OffFileAtomData atoms;
    };
} // namespace cds
#endif
