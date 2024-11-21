#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<size_t>& residueIndices,
                             const std::vector<bool>& residueTER, const PdbFileData& data)
{
    const PdbFileResidueData& residues = data.residues;
    const PdbFileAtomData& atoms       = data.atoms;
    for (size_t n = 0; n < residueIndices.size(); n++)
    {
        size_t residueIndex = residueIndices[n];
        for (size_t atomIndex : residues.atomIndices[residueIndex])
        {
            cds::writeAtomToPdb(stream, residues, residueIndex, atoms, atomIndex);
        }
        if (residueTER[n])
        {
            stream << "TER\n";
        }
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const PdbFileResidueData& residues, size_t residueIndex,
                         const PdbFileAtomData& atoms, size_t atomIndex)
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
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atoms.occupancies[atomIndex];
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atoms.temperatureFactors[atomIndex]
           << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << atoms.elements[atomIndex];
    //    We probably don't want to write charges into pdb file. Width allowed is 2
    stream << std::endl;
}

void cds::writeConectCards(std::ostream& stream, const std::vector<int>& atomNumbers,
                           std::vector<std::pair<size_t, size_t>> connectionIndices)
{ // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but the
  // format is what the format is.
    for (auto& atomPair : connectionIndices)
    {
        stream << "CONECT" << std::right << std::setw(5) << atomNumbers[atomPair.first] << std::right << std::setw(5)
               << atomNumbers[atomPair.second] << "\n";
        stream << "CONECT" << std::right << std::setw(5) << atomNumbers[atomPair.second] << std::right << std::setw(5)
               << atomNumbers[atomPair.first] << "\n";
    }
}
