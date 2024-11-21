#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEDATA_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

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
        std::vector<std::vector<size_t>> atomIndices;
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> chainIds;
        std::vector<std::string> insertionCodes;
    };

    struct PdbFileData
    {
        PdbFileResidueData residues;
        PdbFileAtomData atoms;
    };
} // namespace cds
#endif
