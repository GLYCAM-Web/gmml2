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

cds::AtomPdbData::AtomPdbData(const cds::Atom* atom, std::string recordName_, std::string residueName_,
                              int residueNumber_, std::string chainId_, std::string insertionCode_, double occupancy_,
                              double temperatureFactor_)
    : coordinate({atom->coordinate()}), number({atom->getNumber()}), name({atom->getName()}),
      element({atom->getElement()}), recordName({recordName_}), residueAlternativeLocation({""}),
      residueNumber({residueNumber_}), residueName({residueName_.substr(0, 3)}), chainId({chainId_}),
      insertionCode({insertionCode_}), occupancy({occupancy_}), temperatureFactor({temperatureFactor_})
{}

cds::AtomPdbData::AtomPdbData(std::vector<cds::Atom*> atoms, std::vector<std::string> recordNames,
                              std::vector<std::string> residueNames, std::vector<int> residueNumbers)
    : coordinate(atomCoordinates(atoms)), number(atomNumbers(atoms)), name(atomNames(atoms)),
      element(atomElements(atoms)), recordName(recordNames),
      residueAlternativeLocation(std::vector<std::string>(atoms.size(), "")), residueNumber(residueNumbers),
      residueName(nameSubstrings(residueNames)), chainId(std::vector<std::string>(atoms.size(), "")),
      insertionCode(std::vector<std::string>(atoms.size(), "")), occupancy(std::vector<double>(atoms.size(), 1.0)),
      temperatureFactor(std::vector<double>(atoms.size(), 0.0))
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

cds::ResiduePdbData cds::toResiduePdbData(std::vector<Residue*>& residues)
{
    std::vector<Atom*> atoms;
    std::vector<std::vector<size_t>> indices;
    std::vector<std::string> residueNames;
    std::vector<int> residueNumbers;
    for (auto& residue : residues)
    {
        std::vector<Atom*> residueAtoms = residue->getAtoms();
        indices.push_back(codeUtils::indexVectorWithOffset(atoms.size(), residueAtoms));
        codeUtils::insertInto(atoms, residueAtoms);
        codeUtils::insertInto(residueNames, std::vector<std::string>(residueAtoms.size(), residue->getName()));
        codeUtils::insertInto(residueNumbers, std::vector<int>(residueAtoms.size(), residue->getNumber()));
    }
    std::vector<std::string> recordNames(atoms.size(), "ATOM");
    return ResiduePdbData(indices, AtomPdbData(atoms, recordNames, residueNames, residueNumbers));
}

void cds::writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules)
{
    for (auto& molecule : molecules)
    {
        std::vector<cds::Residue*> residues = molecule->getResidues();
        std::vector<cds::ResidueType> types = residueTypes(residues);
        std::vector<bool> ter               = residueTER(types);
        ResiduePdbData data                 = toResiduePdbData(residues);
        cds::writeMoleculeToPdb(stream, ter, data.atomIndices, data.atomData);
    }
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
            ResiduePdbData data                 = toResiduePdbData(residues);
            cds::writeMoleculeToPdb(stream, ter, data.atomIndices, data.atomData);
        }
        stream << "ENDMDL\n";
    }
}

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<bool>& residueTER,
                             const std::vector<std::vector<size_t>>& indices, const AtomPdbData& atomData)
{
    for (size_t n = 0; n < indices.size(); n++)
    {
        for (size_t index : indices[n])
        {
            cds::writeAtomToPdb(stream, atomData, index);
        }
        if (residueTER[n])
        {
            stream << "TER\n";
        }
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const AtomPdbData& data, size_t index)
{
    const cds::Coordinate& coord = data.coordinate[index];
    stream << std::left << std::setw(6) << data.recordName[index];
    stream << std::right << std::setw(5) << data.number[index] << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << data.name[index];
    stream << std::left << std::setw(1) << data.residueAlternativeLocation[index];
    stream << std::right << std::setw(3) << data.residueName[index] << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << data.chainId[index];
    stream << std::right << std::setw(4) << data.residueNumber[index];
    stream << std::left << std::setw(1) << data.insertionCode[index] << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << data.occupancy[index];
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << data.temperatureFactor[index]
           << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << data.element[index];
    //    We probably don't want to write charges into pdb file. Width allowed is 2
    stream << std::endl;
}

void cds::writeConectCards(std::ostream& stream, std::vector<std::pair<int, int>> connectedAtomNumbers)
{ // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but the
  // format is what the format is.
    for (auto& atomPair : connectedAtomNumbers)
    {
        stream << "CONECT" << std::right << std::setw(5) << atomPair.first << std::right << std::setw(5)
               << atomPair.second << "\n";
        stream << "CONECT" << std::right << std::setw(5) << atomPair.second << std::right << std::setw(5)
               << atomPair.first << "\n";
    }
}
