#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <fstream>
#include <istream>
#include <string>
#include <array>
#include <vector>

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
        return {codeUtils::withoutQuotes(name), codeUtils::withoutQuotes(type), charge, atomIndex};
    }

    void readResidue(std::stringstream& residueStream, lib::LibraryData& data, size_t residueIndex,
                     const std::string& name)
    {
        lib::AtomData& atoms = data.atoms;
        ;
        size_t offset       = atoms.names.size();
        bool hasCoordinates = false;
        std::string line;
        while (getline(residueStream, line) &&
               line.front() !=
                   '!') // Iterate until to the next section that indicates by ! at the beginning of the read line
        {               // Process atom section
            Atom atom = readAtom(line);
            atoms.residueIndex.push_back(residueIndex);
            atoms.names.push_back(atom.name);
            atoms.types.push_back(atom.type);
            atoms.charges.push_back(atom.charge);
            atoms.numbers.push_back(atom.number);
        }
        atoms.coordinates.resize(atoms.names.size(), {0.0, 0.0, 0.0});
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
                    data.bonds.push_back({offset + from - 1, offset + to - 1});
                }
            }
            if (line.find("positions") != std::string::npos)
            { // order of atoms in atom vector corresponds to order in this section
                hasCoordinates = true;
                for (size_t n = offset; n < atoms.names.size(); n++)
                {
                    getline(residueStream, line);
                    std::stringstream ss(line);
                    double x, y, z;
                    ss >> x >> y >> z;
                    atoms.coordinates[n] = {x, y, z};
                }
            }
        }
        data.residues.names.push_back(name);
        data.residues.hasCoordinates.push_back(hasCoordinates);
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
    void parseMolecule(cds::Molecule* molecule, const std::string& filename)
    {
        codeUtils::ensureFileExists(filename);
        std::ifstream fileStream(filename);
        if (fileStream.fail())
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + filename);
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
                codeUtils::RemoveQuotes(line);
                codeUtils::RemoveSpaces(line);
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
        LibraryData data;
        size_t residueCount = residueNames.size();
        std::vector<cds::Residue*> residues;
        residues.reserve(residueCount);
        // Iterate on residue names
        for (size_t residueIndex = 0; residueIndex < residueCount; residueIndex++)
        { // Process the atom section of the file for the corresponding residue
            const std::string& residueName = residueNames[residueIndex];
            std::stringstream residueStream;
            residueStream = extractUnitSection(fileStream, residueName);
            readResidue(residueStream, data, residueIndex, residueName);
            cds::Residue* residue = molecule->addResidue(std::make_unique<cds::Residue>());
            residues.push_back(residue);
            const std::string& name = data.residues.names[residueIndex];
            residue->setName(name);
            residue->determineType(name);
        }
        size_t atomCount = data.atoms.names.size();
        std::vector<cds::Atom*> atoms;
        atoms.reserve(atomCount);
        for (size_t n = 0; n < atomCount; n++)
        {
            size_t residueIndex   = data.atoms.residueIndex[n];
            cds::Residue* residue = residues[residueIndex];
            cds::Atom* atom       = residue->addAtom(std::make_unique<cds::Atom>());
            atoms.push_back(atom);
            atom->setName(data.atoms.names[n]);
            atom->setType(data.atoms.types[n]);
            atom->setCharge(data.atoms.charges[n]);
            atom->setNumber(data.atoms.numbers[n]);
            if (data.residues.hasCoordinates[residueIndex])
            {
                atom->setCoordinate(data.atoms.coordinates[n]);
            }
        }
        for (auto& bond : data.bonds)
        {
            cds::addBond(atoms[bond[0]], atoms[bond[1]]);
        }
    }

} // namespace lib
