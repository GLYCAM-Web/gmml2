#include "include/readers/Pdb/pdbAtom.hpp"

#include "include/readers/Pdb/pdbFunctions.hpp"
#include "include/util/constants.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <string>

namespace gmml
{
    namespace pdb
    {
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        AtomEntry readAtom(const std::string& line)
        {
            std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
            int shift = checkShiftFromSerialNumberOverrun(line);
            int secondShift = checkSecondShiftFromResidueNumberOverrun(line, shift);
            shift += std::max(1, secondShift) - 1;
            double occupancy = 1.0;
            std::string occupancyStr = line.substr(54 + shift, 6);
            try
            {
                occupancy = std::stod(util::RemoveWhiteSpace(occupancyStr));
            }
            catch (...)
            {
                util::log(__LINE__, __FILE__, util::WAR, "Problem converting to occupancy from: " + occupancyStr);
                util::log(__LINE__, __FILE__, util::WAR, "Problematic line is:" + line);
            }
            double temperatureFactor = 0.0;
            std::string temperatureFactorStr = line.substr(60 + shift, 6);
            try
            {
                temperatureFactor = std::stod(util::RemoveWhiteSpace(temperatureFactorStr));
            }
            catch (...)
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::WAR,
                    "Problem converting to temperatureFactor_ from: " + temperatureFactorStr);
                util::log(__LINE__, __FILE__, util::WAR, "Problematic line is:" + line);
            }
            std::string numberStr = line.substr(6, 6 + shift);
            uint number = constants::iNotSet;
            Coordinate coord {0.0, 0.0, 0.0};
            //  In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First
            //  overun for serial doesn't matter as there should be a space between the number and the name. So the
            //  problem is above 999999
            // Dealing with number overruns for serialNumber and residueNumber
            try
            {
                number = std::stoi(util::RemoveWhiteSpace(numberStr));
            }
            catch (...)
            {
                util::log(__LINE__, __FILE__, util::ERR, "Error converting to atom's serial number from: " + numberStr);
            }
            std::string name = util::RemoveWhiteSpace(line.substr(12 + shift, 4));
            name = name.empty() ? "    " : name;
            // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
            if (secondShift > 1)
            {
                shift +=
                    (secondShift -
                     1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
            }
            try
            {
                coord = coordinateFromStrings(
                    util::RemoveWhiteSpace(line.substr(30 + shift, 8)),
                    util::RemoveWhiteSpace(line.substr(38 + shift, 8)),
                    util::RemoveWhiteSpace(line.substr(46 + shift, 8)));
            }
            catch (...)
            {
                util::log(__LINE__, __FILE__, util::ERR, "Error setting coordinate from this line:\n" + line);
            }
            return {recordName, name, number, coord, occupancy, temperatureFactor};
        }
    } // namespace pdb
} // namespace gmml
