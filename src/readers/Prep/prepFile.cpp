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
        void readPrepFile(Molecule* molecule, const std::string& prep_file)
        {
            molecule->setName("prepFile");
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                ReadAllResidues(molecule, in_file);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            SetAtomConnectivities(molecule);
            Generate3dStructures(molecule);
        }

        void readPrepFile(Molecule* molecule, const std::string& prep_file, const std::vector<std::string> queryNames)
        {
            molecule->setName("prepFile");
            util::ensureFileExists(prep_file);
            std::ifstream in_file(prep_file.c_str());
            if (in_file.is_open())
            {
                ReadQueryResidues(molecule, in_file, queryNames);
                in_file.close();
            }
            else
            {
                throw std::runtime_error("Prep file exists but couldn't be opened.");
            }
            // Note that I'm assuming that given a list of query names, the user wants
            // to use all these residues, thus this isn't wasteful:
            util::log(__LINE__, __FILE__, util::INF, "Finished reading prep file. Now setting atomic connectivities.");
            SetAtomConnectivities(molecule);
            util::log(
                __LINE__, __FILE__, util::INF, "Finished setting atomic connectivities. Now generating 3D structures.");
            Generate3dStructures(molecule);
            util::log(__LINE__, __FILE__, util::INF, "Finished, returning from PreFile constructor..");
        }

        //////////////////////////////////////////////////////////
        //                         MUTATORS                     //
        //////////////////////////////////////////////////////////
        void SetAtomConnectivities(Molecule* molecule)
        {
            for (auto& residue : molecule->getResidues())
            {
                util::erratic_cast<PrepResidue*>(residue)->SetConnectivities();
            }
            return;
        }

        void Generate3dStructures(Molecule* molecule)
        {
            for (auto& residue : molecule->getResidues())
            {
                util::erratic_cast<PrepResidue*>(residue)->Generate3dStructure();
            }
            return;
        }

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        void ReadAllResidues(Molecule* molecule, std::istream& in_file)
        {
            std::string line = "";
            getline(in_file, line);
            getline(in_file, line); // first two lines are always blank apparently. smh.
            getline(in_file, line); // This should be first line of residue entry. Title
            while (util::Trim(line).find("STOP") == std::string::npos) /// End of file
            {
                molecule->addResidue(std::make_unique<PrepResidue>(in_file, line));
                getline(in_file, line); // This should be first line of next residue entry or STOP.
            }
        }

        // Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query
        // repeats twice.
        void ReadQueryResidues(Molecule* molecule, std::istream& in_file, const std::vector<std::string>& queryNames)
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
                        molecule->addResidue(std::make_unique<PrepResidue>(in_file, line));
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

        void Write(Molecule* molecule, const std::string& prep_file)
        {
            try
            {
                util::writeToFile(prep_file, [&](std::ostream& stream) { Write(molecule, stream); });
            }
            catch (...)
            {
                throw std::runtime_error("PrepFile could not be created for writing");
            }
        }

        void Write(Molecule* molecule, std::ostream& stream)
        {
            stream << "\n"
                   << "\n";
            for (auto& residue : molecule->getResidues())
            {
                util::erratic_cast<PrepResidue*>(residue)->Write(stream);
            }
            stream << "STOP\n";
        }

        std::string Print(Molecule* molecule)
        {
            std::string out;
            for (auto& residue : molecule->getResidues())
            {
                out += "**********************************************************************************\n";
                out += util::erratic_cast<PrepResidue*>(residue)->toString();
            }
            return out;
        }
    } // namespace prep
} // namespace gmml
