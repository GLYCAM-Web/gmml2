#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    void WriteOffFileUnit(const std::vector<size_t>& residueIndices, const OffFileResidueData& residues,
                          const OffFileAtomData& atoms, std::ostream& stream, const std::string& unitName);
    void WriteResiduesIndividuallyToOffFile(std::ostream& stream, const OffFileData& data);
    void WriteResiduesTogetherToOffFile(std::ostream& stream, const OffFileData& data, const std::string& unitName);
} // namespace cds
#endif
