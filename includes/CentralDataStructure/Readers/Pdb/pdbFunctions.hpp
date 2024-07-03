#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP

#include <string>
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"

namespace pdb
{
    int checkShiftFromSerialNumberOverrun(const std::string& line);
    int checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift = 0);
    cds::Coordinate checkShiftsAndExtractCoordinate(const std::string& line);

}; // namespace pdb
#endif
