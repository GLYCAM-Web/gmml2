#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

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

cds::PdbFileAtomData cds::toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames,
                                            std::vector<double> occupancies, std::vector<double> temperatureFactors)
{
    return {atomCoordinates(atoms), atomNumbers(atoms), atomNames(atoms), atomElements(atoms), recordNames, occupancies,
            temperatureFactors};
}

cds::PdbFileAtomData cds::toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames)
{
    return {atomCoordinates(atoms),
            atomNumbers(atoms),
            atomNames(atoms),
            atomElements(atoms),
            recordNames,
            std::vector<double>(atoms.size(), 1.0),
            std::vector<double>(atoms.size(), 0.0)};
}

cds::PdbFileData cds::toPdbFileData(std::vector<Residue*>& residues)
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
    PdbFileResidueData residueData {indices, residueNumbers(residues), truncatedResidueNames(residues), chainIds,
                                    insertionCodes};
    PdbFileFormat format;
    return PdbFileData {format, {}, residueData, toPdbFileAtomData(atoms, recordNames)};
}

void cds::writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*>& molecules)
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
            PdbFileData data                    = toPdbFileData(residues);
            cds::writeMoleculeToPdb(stream, codeUtils::indexVector(residues), ter, data);
        }
        stream << "ENDMDL\n";
    }
}
