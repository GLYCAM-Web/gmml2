#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include <iomanip> // setw

namespace
{
    std::vector<std::string> nameSubstrings(const std::vector<std::string>& names)
    {
        std::vector<std::string> result;
        result.reserve(names.size());
        for (auto& name : names)
        {
            result.push_back(name.substr(0, 3));
        }
        return result;
    }
} // namespace

cds::AtomPdbData::AtomPdbData(std::vector<cds::Atom*>& atoms_, std::vector<std::string> recordNames_,
                              std::vector<double> occupancies_, std::vector<double> temperatureFactors_)
    : atoms(atoms_), coordinates(atomCoordinates(atoms_)), numbers(atomNumbers(atoms_)), names(atomNames(atoms_)),
      elements(atomElements(atoms_)), recordNames(recordNames_), occupancies(occupancies_),
      temperatureFactors(temperatureFactors_)
{}

cds::AtomPdbData::AtomPdbData(std::vector<cds::Atom*>& atoms_, std::vector<std::string> recordNames_)
    : atoms(atoms_), coordinates(atomCoordinates(atoms_)), numbers(atomNumbers(atoms_)), names(atomNames(atoms_)),
      elements(atomElements(atoms_)), recordNames(recordNames_), occupancies(std::vector<double>(atoms.size(), 1.0)),
      temperatureFactors(std::vector<double>(atoms.size(), 0.0))
{}

cds::ResiduePdbData::ResiduePdbData(std::vector<std::vector<size_t>> atomIndices_, std::vector<int> numbers_,
                                    std::vector<std::string> names_, std::vector<std::string> chainIds_,
                                    std::vector<std::string> insertionCodes_)
    : atomIndices(atomIndices_), numbers(numbers_), names(nameSubstrings(names_)), chainIds(chainIds_),
      insertionCodes(insertionCodes_)
{}

std::vector<bool> cds::residueTER(const std::vector<ResidueType>& types)
{
    std::vector<bool> result;
    result.reserve(types.size());
    for (size_t n = 0; n < types.size(); n++)
    {
        const ResidueType& type = types[n];
        size_t next             = n + 1;
        bool isSugar            = type == cds::ResidueType::Undefined || type == cds::ResidueType::Sugar ||
                       type == cds::ResidueType::Derivative || type == cds::ResidueType::Aglycone;
        bool nextIsCapping           = next < types.size() && types[next] == cds::ResidueType::ProteinCappingGroup;
        bool betweenTwoCappingGroups = (type == cds::ResidueType::ProteinCappingGroup) && nextIsCapping;
        bool isLast                  = next == types.size();
        bool isProtein               = type == cds::ResidueType::Protein;
        result.push_back(isSugar || betweenTwoCappingGroups || (isLast && isProtein));
    }
    return result;
}

cds::PdbWriterData cds::toPdbWriterData(std::vector<Residue*>& residues)
{
    std::vector<Atom*> atoms;
    std::vector<std::vector<size_t>> indices;
    for (auto& residue : residues)
    {
        std::vector<Atom*> residueAtoms = residue->getAtoms();
        indices.push_back(codeUtils::indexVectorWithOffset(atoms.size(), residueAtoms));
        codeUtils::insertInto(atoms, residueAtoms);
    }
    std::vector<std::string> recordNames(atoms.size(), "ATOM");
    std::vector<std::string> chainIds(residues.size(), "");
    std::vector<std::string> insertionCodes(residues.size(), "");
    ResiduePdbData residueData(indices, residueNumbers(residues), residueNames(residues), chainIds, insertionCodes);
    AtomPdbData atomData(atoms, recordNames);
    return PdbWriterData {residueData, atomData};
}

void cds::writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules)
{
    size_t modelCount = molecules.at(0)->getAtoms().at(0)->getNumberOfCoordinateSets();
    for (size_t coordinateSet = 0; coordinateSet < modelCount; coordinateSet++)
    {
        stream << "MODEL " << std::right << std::setw(8) << (coordinateSet + 1) << "\n";
        for (auto& molecule : molecules)
        {
            std::vector<cds::Residue*> residues = molecule->getResidues();
            for (auto& atom : molecule->getAtoms())
            {
                atom->setCurrentCoordinate(coordinateSet);
            }
            std::vector<cds::ResidueType> types = residueTypes(residues);
            std::vector<bool> ter               = residueTER(types);
            PdbWriterData data                  = toPdbWriterData(residues);
            cds::writeMoleculeToPdb(stream, codeUtils::indexVector(residues), ter, data);
        }
        stream << "ENDMDL\n";
    }
}

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<size_t>& residueIndices,
                             const std::vector<bool>& residueTER, const PdbWriterData& data)
{
    const ResiduePdbData& residues = data.residues;
    const AtomPdbData& atoms       = data.atoms;
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
void cds::writeAtomToPdb(std::ostream& stream, const ResiduePdbData& residues, size_t residueIndex,
                         const AtomPdbData& atoms, size_t atomIndex)
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
