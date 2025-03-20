#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepResidue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CodeUtils/logging.hpp" // split
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>
#include <algorithm> // count

using prep::PrepFile;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFile::PrepFile(const std::string& prep_file)
{
    this->setName("prepFile");
    codeUtils::ensureFileExists(prep_file);
    std::ifstream in_file(prep_file.c_str());
    if (in_file.is_open())
    {
        this->ReadAllResidues(in_file);
        in_file.close(); /// Close the prep file
    }
    else
    {
        throw std::runtime_error("Prep file exists but couldn't be opened.");
    }
    this->SetAtomConnectivities();
    this->Generate3dStructures();
}

PrepFile::PrepFile(const std::string& prep_file, const std::vector<std::string> queryNames)
{
    this->setName("prepFile");
    codeUtils::ensureFileExists(prep_file);
    std::ifstream in_file(prep_file.c_str());
    if (in_file.is_open())
    {
        this->ReadQueryResidues(in_file, queryNames);
        in_file.close(); /// Close the prep files
    }
    else
    {
        throw std::runtime_error("Prep file exists but couldn't be opened.");
    }
    // Note that I'm assuming that given a list of query names, the user wants
    // to use all these residues, thus this isn't wasteful:
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished reading prep file. Now setting atomic connectivities.");
    this->SetAtomConnectivities();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished setting atomic connectivities. Now generating 3D structures.");
    this->Generate3dStructures();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished, returning from PreFile constructor..");
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
// std::vector<std::string> PrepFile::GetAllResidueNames()
//{
//	std::vector<std::string> residue_names;
//	for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
//		std::string residue_name = (*it).first;
//		residue_names.push_back(residue_name);
//	}
//	return residue_names;
// }

// gmml::ResidueNameMap PrepFile::GetAllResidueNamesMap()
//{
//	gmml::ResidueNameMap residue_names = gmml::ResidueNameMap();
//	for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
//		std::string residue_name = (*it).first;
//		residue_names[residue_name] = residue_name;
//	}
//	return residue_names;
// }

// std::vector<std::string> PrepFile::GetAllAtomNamesOfResidue(std::string residue_name)
//{
//     std::vector<std::string> atom_names_of_residue;
//     ResidueMap residue_map = GetResidues();
//     prep::PrepResidue* prep_file_residue = residue_map[residue_name];
//     std::vector<PrepFileAtom*> atoms = prep_file_residue->GetAtoms();
//     for(std::vector<PrepFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
//     {
//         PrepFileAtom* atom = (*it);
//         atom_names_of_residue.push_back(atom->GetName());
//     }
//     return atom_names_of_residue;
// }

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void PrepFile::SetAtomConnectivities()
{
    for (auto& residue : this->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->SetConnectivities();
    }
    return;
}

void PrepFile::Generate3dStructures()
{
    for (auto& residue : this->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->Generate3dStructure();
    }
    return;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepFile::ReadAllResidues(std::istream& in_file)
{
    std::string line = "";
    getline(in_file, line);
    getline(in_file, line);                                         // first two lines are always blank apparently. smh.
    getline(in_file, line);                                         // This should be first line of residue entry. Title
    while (codeUtils::Trim(line).find("STOP") == std::string::npos) /// End of file
    {
        this->addResidue(std::make_unique<PrepResidue>(in_file, line));
        getline(in_file, line); // This should be first line of next residue entry or STOP.
    }
}

// Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query repeats twice.
void PrepFile::ReadQueryResidues(std::istream& in_file, const std::vector<std::string>& queryNames)
{
    std::string line = "";
    getline(in_file, line);
    getline(in_file, line); // first two lines of the file are always blank apparently. smh.
    getline(in_file, line); // This should be first line of residue entry. Title.
    while (codeUtils::Trim(line).find("STOP") == std::string::npos) // While not at end of file
    {
        std::streampos firstResidueLinePosition = in_file.tellg(); // save correct position to start reading residue
        std::string savedTitle                  = line;
        // Need to move to line with residue name on it.
        getline(in_file, line);                                                 // blank line
        getline(in_file, line);                                                 // residue name appears here
        std::vector<std::string> residueNameLine = codeUtils::split(line, ' '); // front() string will be name
        if (codeUtils::contains(queryNames, residueNameLine.front()))
        {
            int numberOfTimesToReadInResidue =
                std::count(queryNames.begin(), queryNames.end(), residueNameLine.front());
            //                		std::cout << residueNameLine.front() << " will be read in " <<
            //                numberOfTimesToReadInResidue << " times.\n";
            while (numberOfTimesToReadInResidue > 0)
            {
                in_file.seekg(firstResidueLinePosition); // go back here so the residue constructor works
                line = savedTitle;
                this->addResidue(std::make_unique<PrepResidue>(in_file, line));
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
        // std::cout << "Back out and line is: " << line << std::endl;
    }
    // std::cout << "Ok this is done with line as:\n " << line << std::endl;
}

void PrepFile::Write(const std::string& prep_file)
{
    try
    {
        codeUtils::writeToFile(prep_file,
                               [&](std::ostream& stream)
                               {
                                   this->Write(stream);
                               });
    }
    catch (...)
    {
        throw std::runtime_error("PrepFile could not be created for writing");
    }
}

void PrepFile::Write(std::ostream& stream)
{
    stream << "\n"
           << "\n";
    for (auto& residue : this->getResidues())
    {
        codeUtils::erratic_cast<PrepResidue*>(residue)->Write(stream);
    }
    stream << "STOP\n";
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
std::string PrepFile::Print() const
{
    std::string out;
    for (auto& residue : this->getResidues())
    {
        out += "**********************************************************************************\n";
        out += codeUtils::erratic_cast<PrepResidue*>(residue)->toString();
    }
    return out;
}
