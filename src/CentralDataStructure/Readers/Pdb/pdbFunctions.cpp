#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <iostream>
#include <iomanip>

int pdb::checkShiftFromSerialNumberOverrun(const std::string& line)
{
    int shift = 0;
    if (isdigit(line[11]) && line[20] != ' ')
    {
        shift =
            codeUtils::GetSizeOfIntInString(line.substr(12)); // The shift starts when this gets overrun, not 11. Right?
        std::stringstream ss;
        ss << "Shift of size " << shift << " detected as position 12 is a digit: " << line[11] << " and position 21 >>>"
           << line[20] << "<<< isn't blank: " << line[20] << ":\n           v        v\n"
           << line << "\n";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
    }
    return shift;
}

int pdb::checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift)
{
    return codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
}

cds::Coordinate pdb::checkShiftsAndExtractCoordinate(const std::string& line)
{
    int shift       = pdb::checkShiftFromSerialNumberOverrun(line);
    int secondShift = pdb::checkSecondShiftFromResidueNumberOverrun(line, shift);
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    { // Combine the shifts, but ignore the first shift in residue sequence number.
        shift += (secondShift - 1);
    }
    return pdb::coordinateFromStrings(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)),
                                      codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)),
                                      codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8)));
}

cds::Coordinate pdb::coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z)
{
    try
    {
        return {std::stod(x), std::stod(y), std::stod(z)};
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Could not convert these strings to doubles: " + x + ", " + y + ", " + z + ", ");
        throw;
    }
}

void pdb::expandLine(std::string& line, int length)
{
    int l = line.length();
    if (l < length)
    {
        int space = length - l;
        std::stringstream ss;
        ss << line << std::setw(space) << " ";
        line = ss.str();
    }
}
