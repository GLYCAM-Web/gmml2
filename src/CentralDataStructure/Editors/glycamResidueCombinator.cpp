#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

#include <algorithm>
#include <functional>
#include <cctype>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>

namespace
{
    struct residueMetadata
    {
        std::string isomer          = "";
        std::string resname         = "";
        std::string ringType        = "";
        std::string residueModifier = "";
        std::string configuration   = "";
    };

    residueMetadata Glycam06PrepNameToDetails(const std::string& prepName)
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
        std::string secondLastLetter =
            prepName.substr(prepName.length() - 2, prepName.length() - 2); // Case determines L/D
        std::string glycamCode = secondLastLetter;
        std::transform(glycamCode.begin(), glycamCode.end(), glycamCode.begin(), ::toupper);
        std::cout << "LastLetter is " << lastLetter << "\n";
        std::cout << "glycamCode is " << glycamCode << "\n";

        residueMetadata output;
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
            output.ringType = "p"; // No way to logic this out as the info is lost in this case. So far the only special
                                   // ones have been p.
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

    void removeOMeMethyl(cds::Residue& queryResidue, cds::Atom* OMeOxygen, const std::string& anomericAtomNumber)
    {
        cds::Atom* ch3 = queryResidue.FindAtom("CH3");
        cds::Atom* h31 = queryResidue.FindAtom("H31");
        cds::Atom* h32 = queryResidue.FindAtom("H32");
        cds::Atom* h33 = queryResidue.FindAtom("H33");
        std::cout << "In here OMeOxygen name: " << OMeOxygen->getName() << "\n";
        if (ch3 == nullptr || h31 == nullptr || h32 == nullptr || h33 == nullptr)
        {
            std::string message = "Found methyl oxygen but not the other appropriately named atoms in residue. Glycam "
                                  "combinations cannot be created. This happened for OMeOxygen " +
                                  OMeOxygen->getName() + " in residue: " + queryResidue.getName();
            throw std::runtime_error(message);
        }
        double chargeBeingDeleted = ch3->getCharge() + h31->getCharge() + h32->getCharge() + h33->getCharge();
        OMeOxygen->setCharge(OMeOxygen->getCharge() + chargeBeingDeleted - 0.194);
        queryResidue.deleteAtom(ch3);
        queryResidue.deleteAtom(h31);
        queryResidue.deleteAtom(h32);
        queryResidue.deleteAtom(h33);
        OMeOxygen->setType("Os");
        OMeOxygen->setName("O" + anomericAtomNumber);
        std::cout << "All done here no issues" << std::endl;
        return;
    }

    void removeHydroxyHydrogen(cds::Residue& queryResidue, const std::string hydrogenNumber)
    {
        cds::Atom* hydrogen = queryResidue.FindAtom("H" + hydrogenNumber + "O");
        if (hydrogen == nullptr)
        { // In ROH and some input prep files it's HO1. Otherwise H1O. Fun.
            hydrogen = queryResidue.FindAtom("HO" + hydrogenNumber);
        }
        cds::Atom* oxygen = queryResidue.FindAtom("O" + hydrogenNumber);
        if (hydrogen == nullptr || oxygen == nullptr)
        {
            std::string message =
                "Cannot find appropriately named atoms in residue. Glycam combinations cannot be created. Oxygen "
                "should be "
                "named e.g. O2 and not 2O. Hydrogen to be substituted should be H2O and not HO2. Both must be present. "
                "This may turn into a fatal issue for the atom numbered: " +
                hydrogenNumber + " in residue: " + queryResidue.getName();
            gmml::log(__LINE__, __FILE__, gmml::WAR, message);
            return;
        }
        oxygen->setCharge(oxygen->getCharge() + hydrogen->getCharge() - 0.194);
        queryResidue.deleteAtom(hydrogen);
        oxygen->setType("Os");
        return;
    }

    std::vector<std::string> selectAllAtomsThatCanBeSubstituted(const cds::Residue& queryResidue)
    {
        std::vector<std::string> foundNames;
        std::string delimiter = "";
        for (auto& atom : queryResidue.getAtoms())
        { // if a hydroxyl with a digit in the second position of atom name. e.g. O2, not OHG, and not O1A.
            if (atom->getType() == "Oh" && atom->getNumberFromName() != 0 &&
                isdigit(atom->getName().back())) // If a hydroxyl like O2
            {
                foundNames.push_back(std::to_string(atom->getNumberFromName()));
                std::cout << delimiter << atom->getName();
                delimiter = ",";
            }
        }
        std::cout << "\n";
        return foundNames;
    }
} // namespace

std::vector<std::vector<std::string>> residueCombinator::getCombinations(const std::vector<std::string>& elements)
{
    std::vector<std::vector<std::string>> combinations;
    std::vector<std::string> currentCombination;
    // Recursive function to generate all combinations
    std::function<void(size_t)> generateCombinations;
    generateCombinations = [&](size_t index)
    {
        // Base case: we have reached the end of the input elements,
        // so we add the current combination to the result
        if (index == elements.size())
        {
            combinations.push_back(currentCombination);
            return;
        }
        // Recursive case: we add the element at the current index to
        // the current combination and generate combinations starting
        // at the next index
        currentCombination.push_back(elements[index]);
        generateCombinations(index + 1);
        currentCombination.pop_back();
        // Recursive case: we skip the element at the current index and
        // generate combinations starting at the next index
        generateCombinations(index + 1);
    };
    // Start generating combinations at the first index
    generateCombinations(0);
    combinations.pop_back(); // I don't want the last empty one.
    return combinations;
}

void generateResidueCombination(std::vector<cds::Residue*>& glycamResidueCombinations,
                                const std::vector<std::string> numberCombination, const cds::Residue& templateResidue)
{
    cds::Residue* newResidue = glycamResidueCombinations.emplace_back(new cds::Residue(templateResidue));
    // cds::Residue* newResidue = glycamResidueCombinations.back();
    std::string delimiter    = "";
    std::stringstream numbersAsString;
    for (auto& atomNumber : numberCombination)
    {
        numbersAsString << delimiter << atomNumber;
        delimiter = ",";
        removeHydroxyHydrogen(*newResidue, atomNumber);
        // Set the tail. All can go to the end of the residue. Don't think order matters.
        Atom* atom = newResidue->FindAtom("O" + atomNumber);
        newResidue->moveAtomToLastPosition(atom);
    }
    std::string residueName = GlycamMetadata::GetGlycam06ResidueLinkageCode(numbersAsString.str());
    if (residueName.empty())
    { // Now we have the need for a new residue nomenclature to kick in.
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "No linkage code found for possible combo: " + numbersAsString.str() + " in residue " +
                      templateResidue.getName());
        glycamResidueCombinations.pop_back(); // This should be rare/never?
    }
    else
    {
        residueName += templateResidue.getName().substr(1);
        newResidue->setName(residueName);
        std::cout << numbersAsString.str() << ": " << newResidue->getName() << "\n";
    }
}

// Problem: The input residue will either contain an anomeric oxygen or not.
// Must generate the version that's missing as both are required.
// Combinations that include the anomeric position should include an anomeric oxygen (eg 1GA).
// Combinations that do not include the anomeric position should not (eg 0GA).

void residueCombinator::generateResidueCombinations(std::vector<cds::Residue*>& glycamResidueCombinations,
                                                    const cds::Residue* starterResidue)
{
    // First generate both versions of the residue; with and without anomeric oxygen.
    // One of these gets edited, depending on what's passed in (we don't know yet, need to figure it out in the next
    // steps) Yes this should be two functions. ToDo.
    std::cout << "Copying the input residue " << std::endl;
    cds::Residue residueWithoutAnomericOxygen = *starterResidue;
    cds::Residue residueWithAnomericOxygen    = *starterResidue;
    std::cout << "Guessing anomeric oxygen" << std::endl;
    Atom* anomer = cdsSelections::guessAnomericAtomByInternalNeighbors(residueWithoutAnomericOxygen.getAtoms());
    std::string anomerNumber           = std::to_string(anomer->getNumberFromName());
    std::vector<Atom*> anomerNeighbors = anomer->getNeighbors();
    std::cout << "Names and types of neighbors:" << std::endl;
    for (auto& neighbor : anomerNeighbors)
    {
        std::cout << neighbor->getName() << "_" << neighbor->getType() << std::endl;
    }
    // Lamda functions for the  std::find;
    auto isTypeHydroxy = [](Atom*& a)
    {
        return (a->getType() == "Oh");
    };
    auto isTypeOMe = [](Atom*& a)
    {
        return (a->getType() == "Os" && a->getName() == "O");
    };
    const auto anomericOxygen = std::find_if(anomerNeighbors.begin(), anomerNeighbors.end(), isTypeHydroxy);
    const auto OMeOxygen      = std::find_if(anomerNeighbors.begin(), anomerNeighbors.end(), isTypeOMe);
    if (anomericOxygen != anomerNeighbors.end())
    {
        std::cout << "Anomeric Oxygen Found\n";
        removeHydroxyHydrogen(residueWithoutAnomericOxygen, anomerNumber);
        std::cout << "Copying" << std::endl;
        residueWithAnomericOxygen = residueWithoutAnomericOxygen;
        std::cout << "Deleting" << std::endl;
        residueWithoutAnomericOxygen.deleteAtom(*anomericOxygen);
        std::cout << "Continuing" << std::endl;
    }
    else if (OMeOxygen != anomerNeighbors.end())
    {
        std::cout << "OMe Oxygen Found\n";
        removeOMeMethyl(residueWithoutAnomericOxygen, (*OMeOxygen), anomerNumber);
        std::cout << "Copying" << std::endl;
        residueWithAnomericOxygen = residueWithoutAnomericOxygen;
        std::cout << "Deleting " << (*OMeOxygen)->getName() << std::endl;
        residueWithoutAnomericOxygen.deleteAtom(*OMeOxygen);
        std::cout << "Continuing" << std::endl;
    }
    else
    { // Ok then grow the anomeric oxygen for the residueWithAnomericOxygen
        std::cout << "No Anomeric Oxygen Found in Template\n";
        anomer = cdsSelections::guessAnomericAtomByInternalNeighbors(residueWithAnomericOxygen.getAtoms());
        Coordinate newOxygenCoordinate = cds::coordinateOppositeToNeighborAverage(
            *anomer->getCoordinate(), atomCoordinates(anomer->getNeighbors()), 1.4);
        Atom* newAnomericOxygen = residueWithAnomericOxygen.addAtomToFront(
            std::make_unique<cds::Atom>("O" + anomerNumber, newOxygenCoordinate));
        newAnomericOxygen->setCharge(-0.388);
        newAnomericOxygen->setType("Os");
        addBond(newAnomericOxygen, anomer);
    }
    // Find all positions that can be substituted, ignore the anomer.
    std::vector<std::string> atomNumbers = selectAllAtomsThatCanBeSubstituted(residueWithoutAnomericOxygen);
    // Extra stuff for getting Residue metadata for website that can be turned off:
    // std::ofstream outFileStream;
    // std::string fileName = "latest.pdb";
    // outFileStream.open(fileName.c_str());
    // cds::writeResidueToPdb(outFileStream, starterResidue);
    auto residueInfoToString             = [](const residueMetadata& a)
    {
        return a.isomer + "_" + a.resname + "_" + a.ringType + "_" + a.residueModifier + "_" + a.configuration;
    };
    residueMetadata residueInfo = Glycam06PrepNameToDetails(starterResidue->getName());
    std::cout << "Found name is: " << residueInfoToString(residueInfo) << "\n";
    std::cout << "Anomer: " << anomerNumber << "\n";
    for (auto& atomNumber : atomNumbers)
    {
        std::cout << atomNumber << ",";
    }
    std::cout << "\n";
    // End Extra Stuff.
    // Create the combinations of the numbers
    std::vector<std::vector<std::string>> numberCombinations = residueCombinator::getCombinations(atomNumbers);
    // Create a residue for each of the combinations, copying and modifying the original.
    for (auto& combination : numberCombinations)
    {
        generateResidueCombination(glycamResidueCombinations, combination, residueWithoutAnomericOxygen);
        // ToDo Activate the below when we can handle combinations with anomeric positions.
        // combination.push_back(anomerNumber);
        // generateResidueCombination(glycamResidueCombinations, combination, residueWithAnomericOxygen);
    }
    // Handle the anomeric position separately. No combinations allowed with this position due to limitations with
    // residue naming system.
    cds::Residue* newResidue = glycamResidueCombinations.emplace_back(new cds::Residue(residueWithAnomericOxygen));
    std::string residueName  = anomerNumber;
    residueName              += starterResidue->getName().substr(1);
    newResidue->setName(residueName);
    std::cout << "Added " << residueName << "\n";
    // Write out the 0.. version:
    newResidue  = glycamResidueCombinations.emplace_back(new cds::Residue(residueWithoutAnomericOxygen));
    residueName = "0";
    residueName += starterResidue->getName().substr(1);
    newResidue->setName(residueName);
    std::cout << "Added " << residueName << "\n";
    return;
}
