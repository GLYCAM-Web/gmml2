#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"

#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepResidue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>

void prep::readPrepFile(cds::Molecule* molecule, const std::string& prep_file)
{
    molecule->setName("prepFile");
    codeUtils::ensureFileExists(prep_file);
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

void prep::readPrepFile(
    cds::Molecule* molecule, const std::string& prep_file, const std::vector<std::string> queryNames)
{
    molecule->setName("prepFile");
    codeUtils::ensureFileExists(prep_file);
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
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished reading prep file. Now setting atomic connectivities.");
    SetAtomConnectivities(molecule);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished setting atomic connectivities. Now generating 3D structures.");
    Generate3dStructures(molecule);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished, returning from PreFile constructor..");
}

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void prep::SetAtomConnectivities(cds::Molecule* molecule)
{
    for (auto& residue : molecule->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->SetConnectivities();
    }
    return;
}

void prep::Generate3dStructures(cds::Molecule* molecule)
{
    for (auto& residue : molecule->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->Generate3dStructure();
    }
    return;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void prep::ReadAllResidues(cds::Molecule* molecule, std::istream& in_file)
{
    std::string line = "";
    getline(in_file, line);
    getline(in_file, line);                                         // first two lines are always blank apparently. smh.
    getline(in_file, line);                                         // This should be first line of residue entry. Title
    while (codeUtils::Trim(line).find("STOP") == std::string::npos) /// End of file
    {
        molecule->addResidue(std::make_unique<PrepResidue>(in_file, line));
        getline(in_file, line); // This should be first line of next residue entry or STOP.
    }
}

// Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query repeats twice.
void prep::ReadQueryResidues(cds::Molecule* molecule, std::istream& in_file, const std::vector<std::string>& queryNames)
{
    std::string line = "";
    getline(in_file, line);
    getline(in_file, line); // first two lines of the file are always blank apparently. smh.
    getline(in_file, line); // This should be first line of residue entry. Title.
    while (codeUtils::Trim(line).find("STOP") == std::string::npos) // While not at end of file
    {
        std::streampos firstResidueLinePosition = in_file.tellg(); // save correct position to start reading residue
        std::string savedTitle = line;
        // Need to move to line with residue name on it.
        getline(in_file, line);                                                 // blank line
        getline(in_file, line);                                                 // residue name appears here
        std::vector<std::string> residueNameLine = codeUtils::split(line, ' '); // front() string will be name
        if (codeUtils::contains(queryNames, residueNameLine.front()))
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
            while (codeUtils::Trim(line).find("DONE") == std::string::npos)
            {
                getline(in_file, line);
            }
        }
        getline(in_file, line); // This should be first line of next residue entry or STOP.
    }
}

void prep::Write(cds::Molecule* molecule, const std::string& prep_file)
{
    try
    {
        codeUtils::writeToFile(prep_file, [&](std::ostream& stream) { Write(molecule, stream); });
    }
    catch (...)
    {
        throw std::runtime_error("PrepFile could not be created for writing");
    }
}

void prep::Write(cds::Molecule* molecule, std::ostream& stream)
{
    stream << "\n"
           << "\n";
    for (auto& residue : molecule->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->Write(stream);
    }
    stream << "STOP\n";
}

std::string prep::Print(cds::Molecule* molecule)
{
    std::string out;
    for (auto& residue : molecule->getResidues())
    {
        out += "**********************************************************************************\n";
        out += codeUtils::erratic_cast<PrepResidue*>(residue)->toString();
    }
    return out;
}
