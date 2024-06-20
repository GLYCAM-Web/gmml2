#include "includes/MolecularMetadata/GLYCAM/glycam06PrepToSugarDetail.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include <stdexcept> // throw
#include <iostream>  // cout
#include <algorithm> // for using transform
#include <cctype>    // for using toupper

GlycamMetadata::residueMetadata GlycamMetadata::Glycam06PrepNameToDetails(const std::string& prepName)
{ // This code is a mess because our scientific logic is a mess.
    // Example
    //  input: 0GA
    //  Output: D Glc p  a
    //  Output Isomer resname ringType residueModifier configuration
    //  Note the 0 in 0GA is ignored by this function, GA gives us all the info.
    if (prepName.size() < 3)
    {
        throw std::runtime_error(
            "PrepName is less than 3 characters long, cannot determine the molecular details from GlycamMetadata");
    }
    // Assume is the common situation e.g. 0LB, if not look for weirdos like 0Tv, 0Ko, or 4uA1
    std::string lastLetter =
        prepName.substr(prepName.length() - 1); // If not D,U,A,B, then case determines a(upper) /b (lower) config.
    std::string secondLastLetter = prepName.substr(prepName.length() - 2, prepName.length() - 2); // Case determines L/D
    std::string glycamCode       = secondLastLetter;
    std::transform(glycamCode.begin(), glycamCode.end(), glycamCode.begin(), ::toupper);
    std::cout << "LastLetter is " << lastLetter << "\n";
    std::cout << "glycamCode is " << glycamCode << "\n";

    GlycamMetadata::residueMetadata output;
    (islower(secondLastLetter.at(0))) ? output.isomer = "L" : output.isomer = "D";
    // Configuration Code and ring type from last letter: A/B/U/D
    if (lastLetter == "D")
    {
        output.ringType      = "f";
        output.configuration = "a";
    }
    else if (lastLetter == "U")
    {
        output.ringType      = "f";
        output.configuration = "b";
    }
    else if (lastLetter == "A")
    {
        output.ringType      = "p";
        output.configuration = "a";
    }
    else if (lastLetter == "B")
    {
        output.ringType      = "p";
        output.configuration = "b";
    }
    else // Now need to check last 2 characters and everything is special:
    {
        glycamCode = prepName.substr(prepName.length() - 2);
        std::transform(glycamCode.begin(), glycamCode.end(), glycamCode.begin(), ::toupper);
        (islower(lastLetter.at(0))) ? output.configuration = "b" : output.configuration = "a";
        output.ringType =
            "p"; // No way to logic this out as the info is lost in this case. So far the only special ones have been p.
    }
    std::cout << "So far for " + prepName + " we think isomer ringtype config is " + output.isomer + "_" +
                     output.ringType + "_" + output.configuration + "\n";
    std::cout << "Finding name for code: " + glycamCode + "\n";
    output.resname = GlycamMetadata::GetNameForCode(glycamCode);
    if (output.resname == "") // Dealing with wierdos like 4uA1 or YNP. .
    {
        glycamCode                  = prepName.substr(prepName.length() - 3); // Try the last 3 characters
        output.resname              = GlycamMetadata::GetNameForCode(glycamCode);
        std::string thirdLastLetter = prepName.substr(prepName.length() - 3, prepName.length() - 2);
        (islower(thirdLastLetter.at(0))) ? output.isomer = "L" : output.isomer = "D";
    }
    if (output.resname == "")
    {
        throw std::runtime_error("Did not find a resname in metadata for this prep residue name: " + prepName);
    }

    return output;
}