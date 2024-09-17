#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CodeUtils/containers.hpp"
#include <iomanip> // setw

cds::AtomPdbData cds::toAtomPdbData(const cds::Atom* atom, std::string recordName, std::string residueName,
                                    int residueNumber, std::string chainId, std::string insertionCode, double occupancy,
                                    double temperatureFactor)
{
    std::string residueAlternativeLocation = "";
    return {atom->coordinate(),
            atom->getNumber(),
            atom->getName(),
            atom->getElement(),
            recordName,
            residueAlternativeLocation,
            residueNumber,
            residueName.substr(0, 3),
            chainId,
            insertionCode,
            occupancy,
            temperatureFactor};
}

cds::AtomPdbData cds::toAtomPdbData(const cds::Atom* atom, std::string recordName, std::string residueName,
                                    int residueNumber)
{
    return toAtomPdbData(atom, recordName, residueName, residueNumber, "", "", 1.0, 0.0);
}

std::vector<bool> cds::residueTER(const std::vector<ResidueType>& types)
{
    std::vector<bool> result;
    result.reserve(types.size());
    for (size_t n = 0; n < types.size(); n++)
    {
        auto& type   = types[n];
        size_t next  = n + 1;
        bool isSugar = type == cds::ResidueType::Undefined || type == cds::ResidueType::Sugar ||
                       type == cds::ResidueType::Derivative || type == cds::ResidueType::Aglycone;
        bool nextIsCapping           = next < types.size() && types[next] == cds::ResidueType::ProteinCappingGroup;
        bool betweenTwoCappingGroups = (type == cds::ResidueType::ProteinCappingGroup) && nextIsCapping;
        bool isLast                  = next == types.size();
        bool isProtein               = type == cds::ResidueType::Protein;
        result.push_back(isSugar || betweenTwoCappingGroups || (isLast && isProtein));
    }
    return result;
}

std::vector<cds::AtomPdbData> cds::residuePdbAtoms(Residue* residue)
{
    auto number = residue->getNumber();
    auto name   = residue->getName();
    std::vector<AtomPdbData> atomData;
    auto atoms = residue->getAtoms();
    atomData.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        atomData.push_back(toAtomPdbData(atom, "ATOM", name, number));
    }
    return atomData;
}

std::vector<std::vector<cds::AtomPdbData>> cds::residuePdbAtoms(std::vector<Residue*>& residues)
{
    std::vector<std::vector<AtomPdbData>> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residuePdbAtoms(residue));
    }
    return result;
}

void cds::writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules)
{
    for (auto& molecule : molecules)
    {
        auto residues = molecule->getResidues();
        auto types    = residueTypes(residues);
        auto ter      = residueTER(types);
        auto atoms    = residuePdbAtoms(residues);
        cds::writeMoleculeToPdb(stream, ter, atoms);
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
            auto residues = molecule->getResidues();
            for (auto& atom : molecule->getAtoms())
            {
                atom->setCurrentCoordinate(coordinateSet);
            }
            auto types = residueTypes(residues);
            auto ter   = residueTER(types);
            auto atoms = residuePdbAtoms(residues);
            cds::writeMoleculeToPdb(stream, ter, atoms);
        }
        stream << "ENDMDL\n";
    }
}

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<bool>& residueTER,
                             const std::vector<std::vector<AtomPdbData>>& residueAtoms)
{
    for (size_t n = 0; n < residueAtoms.size(); n++)
    {
        cds::writeResidueToPdb(stream, residueAtoms[n]);
        if (residueTER[n])
        {
            stream << "TER\n";
        }
    }
}

void cds::writeResidueToPdb(std::ostream& stream, const std::vector<AtomPdbData>& atoms)
{
    for (auto& atom : atoms)
    {
        cds::writeAtomToPdb(stream, atom);
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const AtomPdbData& data)
{
    auto coord = data.coordinate;
    stream << std::left << std::setw(6) << data.recordName;
    stream << std::right << std::setw(5) << data.number << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << data.name;
    stream << std::left << std::setw(1) << data.residueAlternativeLocation;
    stream << std::right << std::setw(3) << data.residueName << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << data.chainId;
    stream << std::right << std::setw(4) << data.residueNumber;
    stream << std::left << std::setw(1) << data.insertionCode << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << coord.GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << data.occupancy;
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << data.temperatureFactor << std::left
           << std::setw(10) << " ";
    stream << std::right << std::setw(2) << data.element;
    //    We probably don't want to write charges into pdb file. Width allowed is 2
    stream << std::endl;
    return;
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
