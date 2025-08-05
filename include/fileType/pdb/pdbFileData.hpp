#ifndef INCLUDE_FILETYPE_PDB_PDBFILEDATA_HPP
#define INCLUDE_FILETYPE_PDB_PDBFILEDATA_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/util/formatting.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        struct PdbFileAtomData
        {
            std::vector<Coordinate> coordinates;
            std::vector<uint> numbers;
            std::vector<std::string> names;
            std::vector<std::string> elements;
            std::vector<std::string> recordNames;
            std::vector<double> occupancies;
            std::vector<double> temperatureFactors;
        };

        struct PdbFileResidueData
        {
            std::vector<uint> numbers;
            std::vector<std::string> names;
            std::vector<std::string> chainIds;
            std::vector<std::string> insertionCodes;
        };

        struct PdbFileFormat
        {
            util::floatFormat coordinate = {
                util::textAlignment::right, 8, {3, 3}
            };
            util::floatFormat occupancy = {
                util::textAlignment::right, 6, {2, 2}
            };
            util::floatFormat temperatureFactor = {
                util::textAlignment::right, 6, {2, 2}
            };
        };

        struct PdbFileData
        {
            PdbFileFormat format;
            std::vector<std::string> headerLines;
            PdbFileResidueData residues;
            PdbFileAtomData atoms;
        };
    } // namespace pdb
} // namespace gmml

#endif
