#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    struct AtomOffData
    {
        AtomOffData(std::vector<int> numbers_, std::vector<std::string> names_, std::vector<std::string> types_,
                    std::vector<int> atomicNumbers_, std::vector<double> charges_, std::vector<Coordinate> coordinates_,
                    std::vector<std::vector<int>> children_);
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<int> atomicNumbers;
        std::vector<double> charges;
        std::vector<Coordinate> coordinates;
        std::vector<std::vector<int>> children;
    };

    struct ResidueOffData
    {
        ResidueOffData(std::vector<int> numbers_, std::vector<std::string> names_, std::vector<ResidueType> types_,
                       std::vector<std::vector<size_t>> atomIndices_, std::vector<std::vector<int>> connections_);
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<ResidueType> types;
        std::vector<std::vector<size_t>> atomIndices;
        std::vector<std::vector<int>> atomsConnectedToOtherResidues;
    };

    struct OffWriterData
    {
        ResidueOffData residues;
        AtomOffData atoms;
    };

    OffWriterData toOffWriterData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<cds::Residue*>& residues);
    std::string getOffType(const cds::ResidueType queryType);
    void WriteOffFileUnit(const std::vector<size_t>& residueIndices, const ResidueOffData& residues,
                          const AtomOffData& atoms, std::ostream& stream, const std::string unitName);
    void WriteResiduesIndividuallyToOffFile(std::ostream& stream, const OffWriterData& data);
    void WriteResiduesTogetherToOffFile(std::ostream& stream, const OffWriterData& data, const std::string unitName);
} // namespace cds
#endif
