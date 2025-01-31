#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<cds::Residue*>& residues);
    void WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& indices);
} // namespace cds
#endif
