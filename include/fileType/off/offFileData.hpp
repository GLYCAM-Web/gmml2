#ifndef INCLUDE_FILETYPE_OFF_OFFFILEDATA_HPP
#define INCLUDE_FILETYPE_OFF_OFFFILEDATA_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/residueTypes.hpp"
#include "include/util/formatting.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace off
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
            util::floatFormat coordinate = {
                util::textAlignment::left, 7, {6, 6}
            };
            util::floatFormat charge = {
                util::textAlignment::left, 7, {6, 6}
            };
        };

        struct OffFileData
        {
            OffFileFormat format;
            OffFileResidueData residues;
            OffFileAtomData atoms;
        };
    } // namespace off
} // namespace gmml

#endif
