#include "include/fileType/lib/libraryFile.hpp"

#include "include/util/containers.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <array>
#include <fstream>
#include <istream>
#include <string>
#include <vector>

namespace gmml
{
    namespace
    {
        struct Atom
        {
            std::string name;
            std::string type;
            double charge;
            int number;
        };

        Atom readAtom(const std::string& line)
        {
            int throwaway;
            std::string name;
            std::string type;
            double charge;
            int residueIndex;
            int atomIndex;
            int atomicNumber;
            std::stringstream ss(line);
            ss >> name >> type >> throwaway >> residueIndex >> throwaway >> atomIndex >> atomicNumber >>
                charge; /// Split the line by space to extract attributes of the atom
            return {util::withoutQuotes(name), util::withoutQuotes(type), charge, atomIndex};
        }

        Element toElement(const std::string& name)
        {
            if (!name.empty())
            {
                if (isalpha(name.at(0))) // if first char is in the alphabet
                {
                    return gmml::toElement(name.substr(0, 1)); // return first character as string
                }
            }
            util::log(__LINE__, __FILE__, util::WAR, "Did not find an element for atom named: " + name);
            return Element::Unknown;
        }

        lib::ResidueData readResidue(std::stringstream& residueStream)
        {
            lib::AtomData atoms;
            std::vector<std::array<size_t, 2>> bonds;
            bool hasCoordinates = false;
            std::string line;
            while (getline(residueStream, line) &&
                   line.front() !=
                       '!') // Iterate until to the next section that indicates by ! at the beginning of the read line
            {               // Process atom section
                Atom atom = readAtom(line);
                atoms.names.push_back(atom.name);
                atoms.types.push_back(atom.type);
                atoms.elements.push_back(toElement(atom.name));
                atoms.charges.push_back(atom.charge);
                atoms.numbers.push_back(atom.number);
                atoms.coordinates.push_back({0, 0, 0});
            }
            // Ok atoms are made, now go through the rest of the stream to get their connectivity and coordinates. Not
            // ideal.
            while (getline(residueStream, line))
            {
                if (line.find("unit.connect ") != std::string::npos)
                {
                    while (getline(residueStream, line) && line.front() != '!')
                    {

                        int head_atom_index_;
                        int tail_atom_index_;
                        /// Process connect section
                        std::stringstream ss1(line);
                        ss1 >> head_atom_index_;
                        getline(residueStream, line);
                        std::stringstream ss2(line);
                        ss2 >> tail_atom_index_;
                    }
                }
                if (line.find("connectivity") != std::string::npos)
                {
                    while (getline(residueStream, line) && line.front() != '!')
                    { // Process connectivity section
                        std::stringstream ss(line);
                        size_t from;
                        size_t to;
                        int throwaway;
                        ss >> from >> to >> throwaway;
                        bonds.push_back({from - 1, to - 1});
                    }
                }
                if (line.find("positions") != std::string::npos)
                { // order of atoms in atom vector corresponds to order in this section
                    hasCoordinates = true;
                    for (size_t n = 0; n < atoms.names.size(); n++)
                    {
                        getline(residueStream, line);
                        std::stringstream ss(line);
                        double x, y, z;
                        ss >> x >> y >> z;
                        atoms.coordinates[n] = {x, y, z};
                    }
                }
            }
            return {hasCoordinates, atoms, bonds};
        }

        std::stringstream extractUnitSection(std::istream& inputFileStream, const std::string unitName)
        {
            std::stringstream extractedSection;
            std::string entryLineStart = "!entry." + unitName + ".unit";
            std::string line;
            while ((std::getline(inputFileStream, line)))
            {
                if (line.find("!entry") != std::string::npos && line.find(entryLineStart) == std::string::npos)
                { // If line is an "entry" line that doesn't match this unit, time to leave.
                    // std::cout << "I have never seen this line before in my life: " << line << std::endl;
                    return extractedSection;
                }
                extractedSection << line << std::endl;
            }
            return extractedSection;
        }

    } // namespace

    namespace lib
    {
        LibraryData loadLibraryData(const std::string& filename)
        {
            util::ensureFileExists(filename);
            std::ifstream fileStream(filename);
            if (fileStream.fail())
            {
                util::log(__LINE__, __FILE__, util::ERR, "Could not open this file: " + filename);
                throw std::runtime_error("PdbFile constructor could not open this file: " + filename);
            }
            std::string line;
            // Skip any blank lines at the beginning of the file
            while (line.empty() || line.front() != '!')
            {
                getline(fileStream, line);
            }
            // Read in the array of residue names
            std::vector<std::string> residueNames;
            if (line.find("index") != std::string::npos)
            {
                while (getline(fileStream, line) && line.front() != '!')
                {
                    util::RemoveQuotes(line);
                    util::RemoveSpaces(line);
                    residueNames.push_back(line);
                }
            }
            if (line.find("atoms") ==
                std::string::npos) // If first line passed in for a residue isn't the atoms table, Freak out.
            {
                std::string message =
                    "Error reading library file, I expected the !entry line of an atoms table, but got this: " + line;
                throw std::runtime_error(message);
            }
            std::vector<ResidueData> result;
            size_t residueCount = residueNames.size();
            // Iterate on residue names
            for (size_t residueIndex = 0; residueIndex < residueCount; residueIndex++)
            { // Process the atom section of the file for the corresponding residue
                const std::string& residueName = residueNames[residueIndex];
                std::stringstream residueStream;
                residueStream = extractUnitSection(fileStream, residueName);
                ResidueData residue = readResidue(residueStream);
                result.push_back(residue);
            }
            return {residueNames, result};
        }

    } // namespace lib
} // namespace gmml
