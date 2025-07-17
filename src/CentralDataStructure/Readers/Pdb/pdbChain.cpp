#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/biology.hpp"

#include <string>
#include <vector>
#include <variant>

void pdb::readChain(PdbData& data, size_t moleculeId, std::stringstream& stream_block)
{
    std::string line;
    while (getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if ((recordName == "ATOM") || (recordName == "HETATM"))
        {
            std::stringstream singleResidueSection = extractSingleResidueFromRecordSection(stream_block, line);
            readResidue(data, moleculeId, singleResidueSection, line);
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "In PdbChain Constructor with record that isn't cool: " + recordName);
            break;
        }
    }
    tagTerminalResidues(data, moleculeId);
}

void pdb::tagTerminalResidues(PdbData& data, size_t moleculeId)
{
    size_t nTer = getNTerminal(data, moleculeId);
    if (nTer < data.indices.residueCount)
    {
        data.residues.isNTerminal[nTer] = true;
        data.objects.residues[nTer]->setNTerminal(true);
    }
    size_t cTer = getCTerminal(data, moleculeId);
    if (cTer < data.indices.residueCount)
    {
        data.residues.isCTerminal[cTer] = true;
        data.objects.residues[cTer]->setCTerminal(true);
    }
}

std::stringstream pdb::extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream singleResidueSection;
    pdb::ResidueId residueId(line);
    pdb::ResidueId initialResidueId = residueId;
    while (residueId == initialResidueId)
    {
        singleResidueSection << line << std::endl;
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        if (!std::getline(pdbFileStream, line))
        {
            break; // // If we hit the end, time to leave.
        }
        residueId = ResidueId(line);
    }
    // Go back to previous line position. E.g. was reading HEADER and found TITLE:
    pdbFileStream.seekg(previousLinePosition);
    return singleResidueSection;
}

// Only makes sense for proteins.
// Assumes vector is populated from N-terminal to C-terminal.
size_t pdb::getNTerminal(const PdbData& data, size_t moleculeId)
{
    std::function<bool(const size_t&)> isProtein = [&](size_t id)
    {
        return codeUtils::contains(biology::proteinResidueNames, data.residues.names[id]);
    };
    std::vector<size_t> proteinResidues = codeUtils::vectorFilter(isProtein, data.molecules.residueOrder[moleculeId]);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return data.indices.residueCount;
    }
    return proteinResidues.front();
}

size_t pdb::getCTerminal(const PdbData& data, size_t moleculeId)
{
    std::function<bool(const size_t&)> isProtein = [&](size_t id)
    {
        return codeUtils::contains(biology::proteinResidueNames, data.residues.names[id]);
    };
    std::vector<size_t> proteinResidues = codeUtils::vectorFilter(isProtein, data.molecules.residueOrder[moleculeId]);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return data.indices.residueCount;
    }
    return proteinResidues.back();
}
