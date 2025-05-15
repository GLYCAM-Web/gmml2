#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEDATA_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/residueTypes.hpp"
#include "includes/CodeUtils/formatting.hpp"

#include <vector>
#include <string>

namespace cds
{
    struct OffFileAtomData
    {
        std::vector<uint> numbers;
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<uint> atomicNumbers;
        std::vector<double> charges;
        std::vector<Coordinate> coordinates;
    };

    struct OffFileResidueData
    {
        std::vector<uint> numbers;
        std::vector<std::string> names;
        std::vector<ResidueType> types;
        std::vector<std::vector<size_t>> atomsConnectedToOtherResidues;
    };

    struct OffFileFormat
    {
        codeUtils::floatFormat coordinate = {
            codeUtils::textAlignment::left, 7, {6, 6}
        };
        codeUtils::floatFormat charge = {
            codeUtils::textAlignment::left, 7, {6, 6}
        };
    };

    struct OffFileData
    {
        OffFileFormat format;
        OffFileResidueData residues;
        OffFileAtomData atoms;
    };
} // namespace cds
#endif
