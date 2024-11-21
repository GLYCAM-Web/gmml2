#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <string>
#include <ostream>

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
        std::vector<std::string> types;
        std::vector<std::vector<size_t>> atomIndices;
        std::vector<std::vector<size_t>> atomsConnectedToOtherResidues;
    };

    struct OffFileData
    {
        OffFileResidueData residues;
        OffFileAtomData atoms;
    };

    std::vector<std::string> residueOffTypes(const std::vector<ResidueType>& residues);
    OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<cds::Residue*>& residues);
    void WriteOffFileUnit(const std::vector<size_t>& residueIndices, const OffFileResidueData& residues,
                          const OffFileAtomData& atoms, std::ostream& stream, const std::string& unitName);
    void WriteResiduesIndividuallyToOffFile(std::ostream& stream, const OffFileData& data);
    void WriteResiduesTogetherToOffFile(std::ostream& stream, const OffFileData& data, const std::string& unitName);
} // namespace cds
#endif
