#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    std::vector<std::string> residueOffTypes(const std::vector<ResidueType>& residues);
    OffFileData toOffFileData(const std::vector<Residue*>& residues);
    void serializeResiduesIndividually(std::vector<cds::Residue*>& residues);
} // namespace cds
#endif
