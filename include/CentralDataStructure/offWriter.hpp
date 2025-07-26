#ifndef INCLUDE_CENTRALDATASTRUCTURE_OFFWRITER_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_OFFWRITER_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/off/offFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    off::OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<Residue*>& residues);
    void WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& data);
} // namespace gmml
#endif
