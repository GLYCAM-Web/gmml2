#include "include/CentralDataStructure/Editors/glycamResidueCombinator.hpp"

#include "include/geometry/measurements.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/metadata/glycam06Functions.hpp"
#include "include/readers/Prep/prepDataTypes.hpp"
#include "include/readers/Prep/prepFunctions.hpp"
#include "include/templateGraph/Node.hpp"
#include "include/templateGraph/TotalCycleDecomposition.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace gmml
{
    namespace
    {
        struct residueMetadata
        {
            std::string isomer = "";
            std::string resname = "";
            std::string ringType = "";
            std::string residueModifier = "";
            std::string configuration = "";
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
                throw std::runtime_error("PrepName is less than 3 characters long, cannot determine the molecular "
                                         "details from GlycamMetadata");
            }
            // Assume is the common situation e.g. 0LB, if not look for weirdos like 0Tv, 0Ko, or 4uA1
            std::string lastLetter = prepName.substr(
                prepName.length() - 1); // If not D,U,A,B, then case determines a(upper) /b (lower) config.
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
                output.ringType = "f";
                output.configuration = "a";
            }
            else if (lastLetter == "U")
            {
                output.ringType = "f";
                output.configuration = "b";
            }
            else if (lastLetter == "A")
            {
                output.ringType = "p";
                output.configuration = "a";
            }
            else if (lastLetter == "B")
            {
                output.ringType = "p";
                output.configuration = "b";
            }
            else // Now need to check last 2 characters and everything is special:
            {
                glycamCode = prepName.substr(prepName.length() - 2);
                std::transform(glycamCode.begin(), glycamCode.end(), glycamCode.begin(), ::toupper);
                (islower(lastLetter.at(0))) ? output.configuration = "b" : output.configuration = "a";
                output.ringType = "p"; // No way to logic this out as the info is lost in this case. So far the only
                                       // special ones have been p.
            }
            std::cout << "So far for " + prepName + " we think isomer ringtype config is " + output.isomer + "_" +
                             output.ringType + "_" + output.configuration + "\n";
            std::cout << "Finding name for code: " + glycamCode + "\n";
            output.resname = metadata::GetNameForCode(glycamCode);
            if (output.resname == "") // Dealing with wierdos like 4uA1 or YNP. .
            {
                glycamCode = prepName.substr(prepName.length() - 3); // Try the last 3 characters
                output.resname = metadata::GetNameForCode(glycamCode);
                std::string thirdLastLetter = prepName.substr(prepName.length() - 3, prepName.length() - 2);
                (islower(thirdLastLetter.at(0))) ? output.isomer = "L" : output.isomer = "D";
            }
            if (output.resname == "")
            {
                throw std::runtime_error("Did not find a resname in metadata for this prep residue name: " + prepName);
            }

            return output;
        }

        void removeOMeMethyl(
            prep::PrepData& data, size_t queryResidue, size_t OMeOxygen, const std::string& anomericAtomNumber)
        {
            std::function<size_t(const std::string&)> findAtom = [&](const std::string& name)
            { return prep::findAtom(data, queryResidue, name); };
            std::vector<size_t> toDelete = util::vectorMap(findAtom, {"CH3", "H31", "H32", "H33"});
            size_t atomCount = data.atomCount;
            std::cout << "In here OMeOxygen name: " << data.atoms.name[OMeOxygen] << "\n";
            if (util::contains(toDelete, atomCount))
            {
                std::string message =
                    "Found methyl oxygen but not the other appropriately named atoms in residue. Glycam "
                    "combinations cannot be created. This happened for OMeOxygen " +
                    data.atoms.name[OMeOxygen] + " in residue: " + data.residues.name[queryResidue];
                throw std::runtime_error(message);
            }
            double chargeBeingDeleted = util::vectorSum(0.0, util::indicesToValues(data.atoms.charge, toDelete));
            data.atoms.charge[OMeOxygen] = data.atoms.charge[OMeOxygen] + chargeBeingDeleted - 0.194;
            for (size_t id : toDelete)
            {
                graph::removeNode(data.atomGraph, id);
            }
            data.atoms.type[OMeOxygen] = "Os";
            data.atoms.name[OMeOxygen] = "O" + anomericAtomNumber;
            std::cout << "All done here no issues" << std::endl;
        }

        void removeHydroxyHydrogen(prep::PrepData& data, size_t queryResidue, const std::string& hydrogenNumber)
        {
            std::function<size_t(const std::string&)> findAtom = [&](const std::string& name)
            { return prep::findAtom(data, queryResidue, name); };
            size_t atomCount = data.atomCount;
            size_t hydrogen = findAtom("H" + hydrogenNumber + "O");
            if (hydrogen == atomCount)
            { // In ROH and some input prep files it's HO1. Otherwise H1O. Fun.
                hydrogen = findAtom("HO" + hydrogenNumber);
            }
            size_t oxygen = findAtom("O" + hydrogenNumber);
            if (hydrogen == atomCount || oxygen == atomCount)
            {
                std::string message =
                    "Cannot find appropriately named atoms in residue. Glycam combinations cannot be created. Oxygen "
                    "should be "
                    "named e.g. O2 and not 2O. Hydrogen to be substituted should be H2O and not HO2. Both must be "
                    "present. "
                    "This may turn into a fatal issue for the atom numbered: " +
                    hydrogenNumber + " in residue: " + data.residues.name[queryResidue];
                util::log(__LINE__, __FILE__, util::ERR, message);
                throw std::runtime_error(message);
                return;
            }
            data.atoms.charge[oxygen] = data.atoms.charge[oxygen] + data.atoms.charge[hydrogen] - 0.194;
            data.atoms.type[oxygen] = "Os";
            graph::removeNode(data.atomGraph, hydrogen);
        }

        uint getAtomNumberFromName(const std::string& name)
        {
            std::string convertableNumber = "";
            for (size_t i = 1; i < name.size(); i++)
            {
                if (isdigit(name[i]))
                {
                    convertableNumber += name[i];
                }
            }
            if (convertableNumber.empty())
            {
                return 0;
            }
            return std::stoi(convertableNumber);
        }

        std::vector<std::string> selectAllAtomsThatCanBeSubstituted(const prep::PrepData& data, size_t queryResidue)
        {
            std::vector<std::string> foundNames;
            std::string delimiter = "";
            for (size_t atom : prep::residueAtoms(data, queryResidue))
            { // if a hydroxyl with a digit in the second position of atom name. e.g. O2, not OHG, and not O1A.
                const std::string& name = data.atoms.name[atom];
                uint number = getAtomNumberFromName(name);
                if (data.atoms.type[atom] == "Oh" && number != 0 && !name.empty() &&
                    isdigit(name.back())) // If a hydroxyl like O2
                {
                    foundNames.push_back(std::to_string(number));
                    std::cout << delimiter << name;
                    delimiter = ",";
                }
            }
            std::cout << "\n";
            return foundNames;
        }

        struct AtomIndex : glygraph::Node<AtomIndex>
        {
            AtomIndex() : Node<AtomIndex>("", {}) {};
            size_t atomIndex = 0;
        };

        std::vector<size_t> findCycleAtoms(const prep::PrepData& data, size_t residue)
        {
            std::vector<size_t> atoms = prep::residueAtoms(data, residue);
            std::vector<AtomIndex> graphAtoms(atoms.size());
            for (size_t n = 0; n < atoms.size(); n++)
            {
                graphAtoms[n].atomIndex = atoms[n];
            }
            std::vector<size_t> localIndex = util::indexVector(data.atomCount);
            util::setIndicesTo(localIndex, atoms, util::indexVector(atoms.size()));

            for (size_t edgeId : prep::residueEdges(data, residue))
            {
                const std::array<size_t, 2>& nodes = data.atomGraph.edgeNodes[edgeId];
                graphAtoms[localIndex[nodes[0]]].addChild("", &graphAtoms[localIndex[nodes[1]]]);
            }
            glygraph::Graph<AtomIndex> atomGraph(&graphAtoms[0]);
            std::vector<std::pair<
                std::unordered_set<glygraph::Node<AtomIndex>*>,
                std::unordered_set<glygraph::Edge<AtomIndex>*>>>
                g1Cycles = cycle_decomp::totalCycleDetect(atomGraph);
            std::vector<size_t> cycleAtoms;
            for (std::pair<
                     std::unordered_set<glygraph::Node<AtomIndex>*>,
                     std::unordered_set<glygraph::Edge<AtomIndex>*>> currCyclePair : g1Cycles)
            {
                for (glygraph::Node<AtomIndex>* currAtom : currCyclePair.first)
                {
                    cycleAtoms.push_back(currAtom->getDerivedClass()->atomIndex);
                }
            }
            return cycleAtoms;
        }

        // For now it's the lowest numbered (e.g. C1 lower than C6) ring atom connected to the ring oxygen.
        size_t guessAnomericAtomByInternalNeighbors(const prep::PrepData& data, size_t residue)
        {
            auto isOxygen = [&](size_t n)
            {
                const std::string& name = data.atoms.name[n];
                return !name.empty() && name[0] == 'O';
            };
            std::vector<size_t> cycleAtoms = findCycleAtoms(data, residue);
            size_t cycleOxygen = data.atomCount;
            for (size_t cycleAtom : cycleAtoms)
            {
                if (isOxygen(cycleAtom))
                {
                    cycleOxygen = cycleAtom;
                }
            }
            if (cycleOxygen == data.atomCount)
            {
                throw std::runtime_error(
                    "Did not find a ring oxygen when trying to guess what the anomeric carbon is. Your "
                    "atom names, are likely, strange.");
            }
            std::vector<size_t> anomerCandidates = prep::atomNeighbors(data, cycleOxygen);
            // So deciding this isn't straightforward. For me use case of prep files this is fine. For reading form the
            // PDB I want a meeting first to figure out what we really need. Using a lambda to sort the candidates so
            // that C1 appears before C2 etc.
            std::sort(
                anomerCandidates.begin(),
                anomerCandidates.end(),
                [&](size_t a, size_t b)
                { return (getAtomNumberFromName(data.atoms.name[a]) < getAtomNumberFromName(data.atoms.name[b])); });
            return anomerCandidates[0];
        }
    } // namespace

    std::vector<std::vector<std::string>> getCombinations(const std::vector<std::string>& elements)
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

    Recombined generateResidueCombination(
        prep::PrepData& data,
        const prep::PrepData& reference,
        const std::vector<std::string> numberCombination,
        const size_t templateResidue)
    {
        size_t newResidue = prep::copyResidue(data, reference, templateResidue);
        // Residue* newResidue = glycamResidueCombinations.back();
        std::string delimiter = "";
        std::stringstream numbersAsString;
        for (auto& atomNumber : numberCombination)
        {
            numbersAsString << delimiter << atomNumber;
            delimiter = ",";
            removeHydroxyHydrogen(data, newResidue, atomNumber);
        }
        std::vector<size_t> residueAtoms = prep::residueAtoms(data, newResidue);
        std::vector<size_t> tail;
        for (auto& atomNumber : numberCombination)
        {
            // Set the tail. All can go to the end of the residue. Don't think order matters.
            size_t atom = prep::findAtom(data, newResidue, "O" + atomNumber);
            tail.push_back(atom);
        }
        std::string residueName = metadata::GetGlycam06ResidueLinkageCode(numbersAsString.str());
        if (residueName.empty())
        { // Now we have the need for a new residue nomenclature to kick in.
            util::log(
                __LINE__,
                __FILE__,
                util::WAR,
                "No linkage code found for possible combo: " + numbersAsString.str() + " in residue " +
                    reference.residues.name[templateResidue]); // This should be rare/never?
            return {false, newResidue, residueAtoms, tail};
        }
        else
        {
            data.residues.name[newResidue] = residueName + reference.residues.name[templateResidue].substr(1);
            std::cout << "Added " << newResidue << ": " << data.residues.name[newResidue] << "\n";
            return {true, newResidue, residueAtoms, tail};
        }
    }

    // Problem: The input residue will either contain an anomeric oxygen or not.
    // Must generate the version that's missing as both are required.
    // Combinations that include the anomeric position should include an anomeric oxygen (eg 1GA).
    // Combinations that do not include the anomeric position should not (eg 0GA).

    std::vector<Recombined> generateResidueCombinations(
        prep::PrepData& data, const prep::PrepData& reference, size_t starterResidue)
    {
        std::vector<Recombined> result;
        // First generate both versions of the residue; with and without anomeric oxygen.
        // One of these gets edited, depending on what's passed in (we don't know yet, need to figure it out in the next
        // steps) Yes this should be two functions. ToDo.
        std::cout << "Copying the input residue " << std::endl;
        size_t residueWithoutAnomericOxygen = copyResidue(data, reference, starterResidue);
        size_t residueWithAnomericOxygen = copyResidue(data, reference, starterResidue);
        std::cout << "Guessing anomeric oxygen" << std::endl;
        size_t anomer = guessAnomericAtomByInternalNeighbors(data, residueWithoutAnomericOxygen);
        std::string anomerNumber = std::to_string(getAtomNumberFromName(data.atoms.name[anomer]));
        std::vector<size_t> anomerNeighbors = prep::atomNeighbors(data, anomer);
        std::cout << "Names and types of neighbors:" << std::endl;
        for (size_t neighbor : anomerNeighbors)
        {
            std::cout << data.atoms.name[neighbor] << "_" << data.atoms.type[neighbor] << std::endl;
        }
        auto isTypeHydroxy = [&](size_t n) { return (data.atoms.type[n] == "Oh"); };
        auto isTypeOMe = [&](size_t n) { return (data.atoms.type[n] == "Os" && data.atoms.name[n] == "O"); };
        const auto anomericOxygen = std::find_if(anomerNeighbors.begin(), anomerNeighbors.end(), isTypeHydroxy);
        const auto OMeOxygen = std::find_if(anomerNeighbors.begin(), anomerNeighbors.end(), isTypeOMe);
        if (anomericOxygen != anomerNeighbors.end())
        {
            std::cout << "Anomeric Oxygen Found\n";
            removeHydroxyHydrogen(data, residueWithoutAnomericOxygen, anomerNumber);
            std::cout << "Copying" << std::endl;
            residueWithAnomericOxygen = residueWithoutAnomericOxygen;
            std::cout << "Deleting" << std::endl;
            graph::removeNode(data.atomGraph, *anomericOxygen);
            std::cout << "Continuing" << std::endl;
        }
        else if (OMeOxygen != anomerNeighbors.end())
        {
            std::cout << "OMe Oxygen Found\n";
            removeOMeMethyl(data, residueWithoutAnomericOxygen, *OMeOxygen, anomerNumber);
            std::cout << "Copying" << std::endl;
            residueWithAnomericOxygen = residueWithoutAnomericOxygen;
            std::cout << "Deleting " << data.atoms.name[*OMeOxygen] << std::endl;
            graph::removeNode(data.atomGraph, *OMeOxygen);
            std::cout << "Continuing" << std::endl;
        }
        else
        { // Ok then grow the anomeric oxygen for the residueWithAnomericOxygen
            std::cout << "No Anomeric Oxygen Found in Template\n";
            anomer = guessAnomericAtomByInternalNeighbors(data, residueWithAnomericOxygen);
            Coordinate newOxygenCoordinate = coordinateOppositeToNeighborAverage(
                data.atoms.coordinate[anomer],
                util::indicesToValues(data.atoms.coordinate, prep::atomNeighbors(data, anomer)),
                1.4);

            size_t newAnomericOxygen = prep::addAtom(data, residueWithAnomericOxygen);
            data.atoms.name[newAnomericOxygen] = "O" + anomerNumber;
            data.atoms.type[newAnomericOxygen] = "Os";
            data.atoms.charge[newAnomericOxygen] = -0.388;
            data.atoms.coordinate[newAnomericOxygen] = newOxygenCoordinate;
            graph::addEdge(data.atomGraph, {newAnomericOxygen, anomer});
        }
        // Find all positions that can be substituted, ignore the anomer.
        std::vector<std::string> atomNumbers = selectAllAtomsThatCanBeSubstituted(data, residueWithoutAnomericOxygen);
        // Extra stuff for getting Residue metadata for website that can be turned off:
        // std::ofstream outFileStream;
        // std::string fileName = "latest.pdb";
        // outFileStream.open(fileName.c_str());
        // writeResidueToPdb(outFileStream, starterResidue);
        auto residueInfoToString = [](const residueMetadata& a)
        { return a.isomer + "_" + a.resname + "_" + a.ringType + "_" + a.residueModifier + "_" + a.configuration; };
        residueMetadata residueInfo = Glycam06PrepNameToDetails(reference.residues.name[starterResidue]);
        std::cout << "Found name is: " << residueInfoToString(residueInfo) << "\n";
        std::cout << "Anomer: " << anomerNumber << "\n";
        // End Extra Stuff.
        // Create the combinations of the numbers
        std::vector<std::vector<std::string>> numberCombinations = getCombinations(atomNumbers);
        // Create a residue for each of the combinations, copying and modifying the original.
        for (auto& combination : numberCombinations)
        {
            result.push_back(generateResidueCombination(data, data, combination, residueWithoutAnomericOxygen));
            // ToDo Activate the below when we can handle combinations with anomeric positions.
            // combination.push_back(anomerNumber);
            // generateResidueCombination(data, combination, residueWithAnomericOxygen);
        }
        // Handle the anomeric position separately. No combinations allowed with this position due to limitations with
        // residue naming system.
        auto addSpecialCase = [&](size_t templateResidue, const std::string& prefix, const std::string& displayStr)
        {
            size_t newResidue = copyResidue(data, data, templateResidue);
            std::string name = prefix + reference.residues.name[starterResidue].substr(1);
            data.residues.name[newResidue] = name;
            result.push_back({true, newResidue, prep::residueAtoms(data, newResidue), {}});
            std::cout << "Added " << displayStr << " " << newResidue << ": " << name << "\n";
        };
        addSpecialCase(residueWithAnomericOxygen, anomerNumber, "with anomeric");
        addSpecialCase(residueWithoutAnomericOxygen, "0", "without anomeric");
        return result;
    }
} // namespace gmml
