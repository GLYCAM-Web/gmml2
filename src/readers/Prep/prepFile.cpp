#include "include/readers/Prep/prepFile.hpp"

#include "include/CentralDataStructure/molecule.hpp"
#include "include/readers/Prep/prepAtom.hpp"
#include "include/readers/Prep/prepResidue.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>

namespace gmml
{
    namespace prep
    {
        void readPrepFile(
            Molecule* molecule, std::vector<PrepResidueProperties>& properties, const std::string& prep_file)
        {
            molecule->setName("prepFile");
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                readAllResidues(molecule, properties, in_file);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            setAtomConnectivities(molecule, properties);
            generate3dStructures(molecule, properties);
        }

        void readPrepFile(
            Molecule* molecule,
            std::vector<PrepResidueProperties>& properties,
            const std::string& prep_file,
            const std::vector<std::string> queryNames)
        {
            molecule->setName("prepFile");
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                readQueryResidues(molecule, properties, in_file, queryNames);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            // Note that I'm assuming that given a list of query names, the user wants
            // to use all these residues, thus this isn't wasteful:
            util::log(__LINE__, __FILE__, util::INF, "Finished reading prep file. Now setting atomic connectivities.");
            setAtomConnectivities(molecule, properties);
            util::log(
                __LINE__, __FILE__, util::INF, "Finished setting atomic connectivities. Now generating 3D structures.");
            generate3dStructures(molecule, properties);
            util::log(__LINE__, __FILE__, util::INF, "Finished, returning from PreFile constructor..");
        }

        //////////////////////////////////////////////////////////
        //                         MUTATORS                     //
        //////////////////////////////////////////////////////////
        void setAtomConnectivities(Molecule* molecule, const std::vector<PrepResidueProperties>& properties)
        {
            std::vector<Residue*> residues = molecule->getResidues();
            for (size_t n = 0; n < residues.size(); n++)
            {
                setConnectivities(residues[n], properties[n]);
            }
        }

        void generate3dStructures(Molecule* molecule, std::vector<PrepResidueProperties>& properties)
        {
            std::vector<Residue*> residues = molecule->getResidues();
            for (size_t n = 0; n < residues.size(); n++)
            {
                generate3dStructure(residues[n], properties[n]);
            }
        }

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        void readAllResidues(Molecule* molecule, std::vector<PrepResidueProperties>& properties, std::istream& in_file)
        {
            std::string line = "";
            getline(in_file, line);
            getline(in_file, line); // first two lines are always blank apparently. smh.
            getline(in_file, line); // This should be first line of residue entry. Title
            while (util::Trim(line).find("STOP") == std::string::npos) /// End of file
            {
                Residue* residue = molecule->addResidue(std::make_unique<Residue>());
                properties.push_back(PrepResidueProperties());
                initializePrepResidue(residue, properties.back(), in_file, line);
                getline(in_file, line); // This should be first line of next residue entry or STOP.
            }
        }

        // Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query
        // repeats twice.
        void readQueryResidues(
            Molecule* molecule,
            std::vector<PrepResidueProperties>& properties,
            std::istream& in_file,
            const std::vector<std::string>& queryNames)
        {
            std::string line = "";
            getline(in_file, line);
            getline(in_file, line); // first two lines of the file are always blank apparently. smh.
            getline(in_file, line); // This should be first line of residue entry. Title.
            while (util::Trim(line).find("STOP") == std::string::npos) // While not at end of file
            {
                std::streampos firstResidueLinePosition =
                    in_file.tellg(); // save correct position to start reading residue
                std::string savedTitle = line;
                // Need to move to line with residue name on it.
                getline(in_file, line);                                            // blank line
                getline(in_file, line);                                            // residue name appears here
                std::vector<std::string> residueNameLine = util::split(line, ' '); // front() string will be name
                if (util::contains(queryNames, residueNameLine.front()))
                {
                    int numberOfTimesToReadInResidue =
                        std::count(queryNames.begin(), queryNames.end(), residueNameLine.front());
                    while (numberOfTimesToReadInResidue > 0)
                    {
                        in_file.seekg(firstResidueLinePosition); // go back here so the residue constructor works
                        line = savedTitle;
                        Residue* residue = molecule->addResidue(std::make_unique<Residue>());
                        properties.push_back(PrepResidueProperties());
                        initializePrepResidue(residue, properties.back(), in_file, line);
                        --numberOfTimesToReadInResidue;
                    }
                }
                else
                { // need to flush the lines until we find a residue we want.
                    while (util::Trim(line).find("DONE") == std::string::npos)
                    {
                        getline(in_file, line);
                    }
                }
                getline(in_file, line); // This should be first line of next residue entry or STOP.
            }
        }

        void write(
            Molecule* molecule, const std::vector<PrepResidueProperties>& properties, const std::string& prep_file)
        {
            try
            {
                util::writeToFile(prep_file, [&](std::ostream& stream) { write(molecule, properties, stream); });
            }
            catch (...)
            {
                throw std::runtime_error("PrepFile could not be created for writing");
            }
        }

        void write(Molecule* molecule, const std::vector<PrepResidueProperties>& properties, std::ostream& stream)
        {
            stream << "\n"
                   << "\n";
            std::vector<Residue*> residues = molecule->getResidues();
            for (size_t n = 0; n < residues.size(); n++)
            {
                write(residues[n], properties[n], stream);
            }
            stream << "STOP\n";
        }

        std::string print(Molecule* molecule, const std::vector<PrepResidueProperties>& properties)
        {
            std::string out;
            std::vector<Residue*> residues = molecule->getResidues();
            for (size_t n = 0; n < residues.size(); n++)
            {
                out += "**********************************************************************************\n";
                out += toString(residues[n], properties[n]);
            }
            return out;
        }
    } // namespace prep
} // namespace gmml
