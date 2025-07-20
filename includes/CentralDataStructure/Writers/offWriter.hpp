#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace cds
{
    OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<cds::Residue*>& residues);
    void WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& data);
} // namespace cds
#endif
