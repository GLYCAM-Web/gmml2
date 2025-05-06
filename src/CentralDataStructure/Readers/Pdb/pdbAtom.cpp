#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Writers/print.hpp"
#include "includes/CodeUtils/constants.hpp" // gmml::iNotSet
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
void pdb::readAtom(cds::Atom* atom, const std::string& line)
{
    //  In the PDB file the residue number overruns after 9999 and serial number overruns after 99999. First overun for
    //  serial doesn't matter as there should be a space between the number and the name. So the problem is above 999999
    // Dealing with number overruns for serialNumber and residueNumber
    int shift = pdb::checkShiftFromSerialNumberOverrun(line);
    try
    {
        atom->setNumber(std::stoi(codeUtils::RemoveWhiteSpace(line.substr(6, 6 + shift))));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Error converting to atom's serial number from: " + line.substr(6, 6 + shift));
        atom->setNumber(constants::iNotSet);
    }
    std::string atomName = codeUtils::RemoveWhiteSpace(line.substr(12 + shift, 4));
    if (atomName.empty())
    {
        atom->setName("    ");
    }
    else
    {
        atom->setName(atomName);
    }
    int secondShift = pdb::checkSecondShiftFromResidueNumberOverrun(line, shift);
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    {
        shift +=
            (secondShift - 1); // Combine the shifts, but ignore the first shift in residue sequence number from here on
    }
    try
    {
        atom->setCoordinate(pdb::coordinateFromStrings(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)),
                                                       codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)),
                                                       codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8))));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error setting coordinate from this line:\n" + line);
    }
}
