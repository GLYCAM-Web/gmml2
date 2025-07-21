#ifndef INCLUDES_WRITERS_OFFWRITER_HPP
#define INCLUDES_WRITERS_OFFWRITER_HPP

#include "include/CentralDataStructure/FileFormats/offFileData.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<Residue*>& residues);
    void WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& data);
} // namespace gmml
#endif
