#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace

#include <string>

pdb::ResidueId pdb::pdbResidueId(const PdbData& data, size_t residueId)
{
    cds::Residue* residue = data.indices.residues[residueId];
    return ResidueId(residue->getName(), std::to_string(residue->getNumber()), data.residues.insertionCodes[residueId],
                     data.residues.chainIds[residueId]);
}

std::string pdb::residueStringId(const PdbData& data, size_t residueId)
{
    return pdbResidueId(data, residueId).print();
}

size_t pdb::addPdbAtom(PdbData& data, size_t residueId, const std::string& line)
{
    int shift                = checkShiftFromSerialNumberOverrun(line);
    int secondShift          = checkSecondShiftFromResidueNumberOverrun(line, shift);
    shift                    += std::max(1, secondShift) - 1;
    double occupancy         = 1.0;
    std::string occupancyStr = line.substr(54 + shift, 6);
    try
    {
        occupancy = std::stod(codeUtils::RemoveWhiteSpace(occupancyStr));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problem converting to occupancy from: " + occupancyStr);
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
    }
    double temperatureFactor         = 0.0;
    std::string temperatureFactorStr = line.substr(60 + shift, 6);
    try
    {
        temperatureFactor = std::stod(codeUtils::RemoveWhiteSpace(temperatureFactorStr));
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Problem converting to temperatureFactor_ from: " + temperatureFactorStr);
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Problematic line is:" + line);
    }
    cds::Atom* atom = data.indices.residues[residueId]->addAtom(std::make_unique<cds::Atom>());
    readAtom(atom, line);
    size_t atomId = data.indices.atoms.size();
    data.indices.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.atoms.recordNames.push_back(codeUtils::RemoveWhiteSpace(line.substr(0, 6)));
    data.atoms.occupancies.push_back(occupancy);
    data.atoms.temperatureFactors.push_back(temperatureFactor);
    return atomId;
}

size_t pdb::addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& c)
{
    cds::Atom* atom = data.indices.residues[residueId]->addAtom(std::make_unique<cds::Atom>(name, c));
    size_t atomId   = data.indices.atoms.size();
    data.indices.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.atoms.recordNames.push_back("ATOM");
    data.atoms.occupancies.push_back(1.0);
    data.atoms.temperatureFactors.push_back(0.0);
    return atomId;
}

void pdb::deletePdbAtom(PdbData& data, size_t residueId, size_t atomId)
{
    if (atomId < data.indices.atoms.size())
    {
        cds::Atom* atom = data.indices.atoms[atomId];
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Deleting atom with id: " + atom->getName() + "_" + std::to_string(atom->getNumber()));
        codeUtils::eraseNth(atomId, data.atoms.recordNames);
        codeUtils::eraseNth(atomId, data.atoms.occupancies);
        codeUtils::eraseNth(atomId, data.atoms.temperatureFactors);
        codeUtils::eraseNth(atomId, data.indices.atomResidue);
        codeUtils::eraseNth(atomId, data.indices.atoms);
        data.indices.residues[residueId]->deleteAtom(atom);
    }
}

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
void pdb::readResidue(PdbData& data, size_t residueId, std::stringstream& singleResidueSecion, std::string firstLine)
{
    cds::Residue* residue = data.indices.residues[residueId];
    ResidueId resId(firstLine);
    residue->setName(resId.getName());
    residue->setNumber(std::stoi(resId.getNumber()));
    data.residues.insertionCodes.push_back(resId.getInsertionCode());
    data.residues.chainIds.push_back(resId.getChainId());
    data.residues.hasTerCard.push_back(false);
    std::string firstFoundAlternativeLocationIndicator = resId.getAlternativeLocation(); // Normally empty
    std::string line;
    while (getline(singleResidueSecion, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if ((recordName == "ATOM") || (recordName == "HETATM"))
        {
            ResidueId id(line); // Check alternativeLocation and ignore any that aren't the same as the first one.
            if (id.getAlternativeLocation().empty() ||
                id.getAlternativeLocation() == firstFoundAlternativeLocationIndicator)
            { // If no alternative location (normal case) or alternateLocation is the first one (normally "A").
                addPdbAtom(data, residueId, line);
            }
            else if (firstFoundAlternativeLocationIndicator.empty())
            { // Sometimes first atom has one location, but later atoms have alternatives. Just set and use the first
              // one (normally "A")
                firstFoundAlternativeLocationIndicator = id.getAlternativeLocation();
                addPdbAtom(data, residueId, line);
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Skipping line with alternative location: " + line);
            }
        }
    }
    residue->SetType(residue->determineType(residue->getName()));
}

void pdb::readResidue(PdbData& data, size_t residueId, size_t referenceResidue)
{ // should instead call copy constructor and then rename with residueName?
    cds::Residue* residue = data.indices.residues[residueId];
    residue->SetType(residue->determineType(residue->getName()));
    data.residues.insertionCodes.push_back(data.residues.insertionCodes[referenceResidue]);
    data.residues.chainIds.push_back(data.residues.chainIds[referenceResidue]);
    data.residues.hasTerCard.push_back(false);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string pdb::getNumberAndInsertionCode(const PdbData& data, size_t residueId)
{
    return std::to_string(data.indices.residues[residueId]->getNumber()) + data.residues.insertionCodes[residueId];
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void pdb::modifyNTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "NH3+")
    {
        size_t atomId = codeUtils::indexOf(data.indices.atoms, data.indices.residues[residueId]->FindAtom("H"));
        deletePdbAtom(data, residueId, atomId);
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void pdb::modifyCTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    auto findAtom = [&](const std::string& name)
    {
        return findResidueAtom(data, residueId, name);
    };
    auto coord = [&](size_t n)
    {
        return data.indices.atoms[n]->coordinate();
    };
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "CO2-")
    {
        size_t atomOXT = findAtom("OXT");
        if (atomOXT == data.indices.atoms.size())
        {
            // I don't like this, but at least it's somewhat contained:
            size_t atomCA            = findAtom("CA");
            size_t atomC             = findAtom("C");
            size_t atomO             = findAtom("O");
            cds::Atom* atomOptr      = data.indices.atoms[atomO];
            cds::Coordinate oxtCoord = cds::calculateCoordinateFromInternalCoords(coord(atomCA), coord(atomC),
                                                                                  coord(atomO), 120.0, 180.0, 1.25);
            addPdbAtom(data, residueId, "OXT", oxtCoord);
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Created new atom named OXT after " + atomOptr->getName() + "_" +
                          std::to_string(atomOptr->getNumber()));
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}
