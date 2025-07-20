#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"

#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
pdb::AtomEntry pdb::readAtom(const std::string& line)
{
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    int shift = checkShiftFromSerialNumberOverrun(line);
    int secondShift = checkSecondShiftFromResidueNumberOverrun(line, shift);
    shift += std::max(1, secondShift) - 1;
    double occupancy = 1.0;
    std::string occupancyStr = line.substr(54 + shift, 6);
    try
    {
        occupancy = std::stod(codeUtils::RemoveWhiteSpace(occupancyStr));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problem converting to occupancy from: " + occupancyStr);
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
    }
    double temperatureFactor = 0.0;
    std::string temperatureFactorStr = line.substr(60 + shift, 6);
    try
    {
        temperatureFactor = std::stod(codeUtils::RemoveWhiteSpace(temperatureFactorStr));
    }
    catch (...)
    {
        gmml::log(
            __LINE__, __FILE__, gmml::WAR, "Problem converting to temperatureFactor_ from: " + temperatureFactorStr);
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
    }
    std::string numberStr = line.substr(6, 6 + shift);
    uint number = constants::iNotSet;
    cds::Coordinate coord {0.0, 0.0, 0.0};
    //  In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First overun for
    //  serial doesn't matter as there should be a space between the number and the name. So the problem is above 999999
    // Dealing with number overruns for serialNumber and residueNumber
    try
    {
        number = std::stoi(codeUtils::RemoveWhiteSpace(numberStr));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error converting to atom's serial number from: " + numberStr);
    }
    std::string name = codeUtils::RemoveWhiteSpace(line.substr(12 + shift, 4));
    name = name.empty() ? "    " : name;
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    {
        shift +=
            (secondShift - 1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
    }
    try
    {
        coord = pdb::coordinateFromStrings(
            codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)),
            codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)),
            codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8)));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting coordinate from this line:\n" + line);
    }
    return {recordName, name, number, coord, occupancy, temperatureFactor};
}
