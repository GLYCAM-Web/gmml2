#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/formatting.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <array>
#include <vector>

void cds::writeAssemblyToPdb(std::ostream& stream, const assembly::Graph& graph,
                             const std::vector<std::vector<size_t>>& residueIndices,
                             const std::vector<std::vector<bool>>& residueTER,
                             const std::vector<std::array<size_t, 2>>& connectionIndices, const PdbFileData& data)
{
    for (auto& line : data.headerLines)
    {
        stream << "HEADER    " << line << "\n";
    }
    for (size_t n = 0; n < residueIndices.size(); n++)
    {
        cds::writeMoleculeToPdb(stream, graph, residueIndices[n], residueTER[n], data);
    }
    cds::writeConectCards(stream, data.atoms.numbers, connectionIndices);
}

void cds::writeMoleculeToPdb(std::ostream& stream, const assembly::Graph& graph,
                             const std::vector<size_t>& residueIndices, const std::vector<bool>& residueTER,
                             const PdbFileData& data)
{
    const PdbFileResidueData& residues = data.residues;
    const PdbFileAtomData& atoms       = data.atoms;
    for (size_t n = 0; n < residueIndices.size(); n++)
    {
        size_t residueIndex = residueIndices[n];
        for (size_t atomIndex : graph.residues.nodes.elements[residueIndex])
        {
            cds::writeAtomToPdb(stream, data.format, residues, residueIndex, atoms, atomIndex);
        }
        if (residueTER[n])
        {
            stream << "TER\n";
        }
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const PdbFileFormat& format, const PdbFileResidueData& residues,
                         size_t residueIndex, const PdbFileAtomData& atoms, size_t atomIndex)
{
    std::string residueAlternativeLocation = ""; // don't know if this should live in residue or atom..
    const cds::Coordinate& coord           = atoms.coordinates[atomIndex];
    stream << std::left << std::setw(6) << atoms.recordNames[atomIndex];
    stream << std::right << std::setw(5) << atoms.numbers[atomIndex] << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << atoms.names[atomIndex];
    stream << std::left << std::setw(1) << residueAlternativeLocation;
    stream << std::right << std::setw(3) << residues.names[residueIndex] << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << residues.chainIds[residueIndex];
    stream << std::right << std::setw(4) << residues.numbers[residueIndex];
    stream << std::left << std::setw(1) << residues.insertionCodes[residueIndex] << std::left << std::setw(3) << " ";
    codeUtils::writeFloat(stream, format.coordinate, coord.GetX());
    codeUtils::writeFloat(stream, format.coordinate, coord.GetY());
    codeUtils::writeFloat(stream, format.coordinate, coord.GetZ());
    codeUtils::writeFloat(stream, format.occupancy, atoms.occupancies[atomIndex]);
    codeUtils::writeFloat(stream, format.temperatureFactor, atoms.temperatureFactors[atomIndex]);
    stream << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << atoms.elements[atomIndex];
    //    We probably don't want to write charges into pdb file. Width allowed is 2
    stream << std::endl;
}

void cds::writeConectCards(std::ostream& stream, const std::vector<int>& atomNumbers,
                           const std::vector<std::array<size_t, 2>>& connectionIndices)
{ // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but the
  // format is what the format is.
    for (auto& atomPair : connectionIndices)
    {
        auto writeLine = [&](int a, int b)
        {
            stream << "CONECT" << std::right << std::setw(5) << a << std::right << std::setw(5) << b << "\n";
        };
        int first  = atomNumbers[atomPair[0]];
        int second = atomNumbers[atomPair[1]];
        writeLine(first, second);
        writeLine(second, first);
    }
}
