#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEDATA_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/formatting.hpp"

#include <string>
#include <vector>

namespace cds
{
    struct PdbFileAtomData
    {
        std::vector<Coordinate> coordinates;
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> elements;
        std::vector<std::string> recordNames;
        std::vector<double> occupancies;
        std::vector<double> temperatureFactors;
    };

    struct PdbFileResidueData
    {
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> chainIds;
        std::vector<std::string> insertionCodes;
    };

    struct PdbFileFormat
    {
        codeUtils::floatFormat coordinate = {
            codeUtils::textAlignment::right, 8, {3, 3}
        };
        codeUtils::floatFormat occupancy = {
            codeUtils::textAlignment::right, 6, {2, 2}
        };
        codeUtils::floatFormat temperatureFactor = {
            codeUtils::textAlignment::right, 6, {2, 2}
        };
    };

    struct PdbFileData
    {
        PdbFileFormat format;
        std::vector<std::string> headerLines;
        PdbFileResidueData residues;
        PdbFileAtomData atoms;
    };
} // namespace cds
#endif
