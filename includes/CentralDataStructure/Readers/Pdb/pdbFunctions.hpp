#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <string>

namespace pdb
{
    int checkShiftFromSerialNumberOverrun(const std::string& line);
    int checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift = 0);
    cds::Coordinate checkShiftsAndExtractCoordinate(const std::string& line);
    cds::Coordinate coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z);
    void expandLine(std::string& line, int length);
    void print(const cds::Coordinate& coord, std::ostream& out);
}; // namespace pdb
#endif
