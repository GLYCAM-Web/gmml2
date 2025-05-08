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

using pdb::PdbResidue;

pdb::ResidueId pdb::pdbResidueId(const PdbData& data, size_t residueId)
{
    PdbResidue* residue = codeUtils::erratic_cast<PdbResidue*>(data.indices.residues[residueId]);
    return ResidueId(residue->getName(), std::to_string(residue->getNumber()), residue->getInsertionCode(),
                     residue->getChainId());
}

cds::Atom* PdbResidue::addPdbAtom(PdbData& data, size_t residueId, const std::string& line)
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
    cds::Atom* atom = this->addAtom(std::make_unique<cds::Atom>());
    readAtom(atom, line);
    data.indices.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.atoms.recordNames.push_back(codeUtils::RemoveWhiteSpace(line.substr(0, 6)));
    data.atoms.occupancies.push_back(occupancy);
    data.atoms.temperatureFactors.push_back(temperatureFactor);
    return atom;
}

cds::Atom* PdbResidue::addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& c)
{
    cds::Atom* atom = this->addAtom(std::make_unique<cds::Atom>(name, c));
    data.indices.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.atoms.recordNames.push_back("ATOM");
    data.atoms.occupancies.push_back(1.0);
    data.atoms.temperatureFactors.push_back(0.0);
    return atom;
}

void PdbResidue::deletePdbAtom(PdbData& data, size_t, cds::Atom* atom)
{
    if (atom != nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Deleting atom with id: " + atom->getName() + "_" + std::to_string(atom->getNumber()));
        size_t index = codeUtils::indexOf(data.indices.atoms, atom);
        codeUtils::eraseNth(index, data.atoms.recordNames);
        codeUtils::eraseNth(index, data.atoms.occupancies);
        codeUtils::eraseNth(index, data.atoms.temperatureFactors);
        codeUtils::eraseNth(index, data.indices.atomResidue);
        codeUtils::eraseNth(index, data.indices.atoms);
        this->deleteAtom(atom);
    }
}

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidue::PdbResidue(PdbData& data, size_t residueId, std::stringstream& singleResidueSecion, std::string firstLine)
{
    this->setPdbResidue();
    ResidueId resId(firstLine);
    this->setName(resId.getName());
    this->setNumber(std::stoi(resId.getNumber()));
    this->setInsertionCode(resId.getInsertionCode());
    this->setChainId(resId.getChainId());
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
                this->addPdbAtom(data, residueId, line);
            }
            else if (firstFoundAlternativeLocationIndicator.empty())
            { // Sometimes first atom has one location, but later atoms have alternatives. Just set and use the first
              // one (normally "A")
                firstFoundAlternativeLocationIndicator = id.getAlternativeLocation();
                this->addPdbAtom(data, residueId, line);
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Skipping line with alternative location: " + line);
            }
        }
    }
    this->SetType(this->determineType(this->getName()));
}

PdbResidue::PdbResidue(PdbData&, size_t, const std::string residueName, const PdbResidue* referenceResidue)
    : cds::Residue(residueName, referenceResidue)
{ // should instead call copy constructor and then rename with residueName?
    this->setPdbResidue();
    this->setInsertionCode(referenceResidue->getInsertionCode());
    this->setChainId(referenceResidue->getChainId());
    this->SetType(this->determineType(this->getName()));
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
const std::string PdbResidue::getNumberAndInsertionCode() const
{
    return std::to_string(this->getNumber()) + this->getInsertionCode();
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void PdbResidue::modifyNTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "NH3+")
    {
        this->deletePdbAtom(data, residueId, this->FindAtom("H"));
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void PdbResidue::modifyCTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "CO2-")
    {
        const cds::Atom* atom = this->FindAtom("OXT");
        if (atom == nullptr)
        {
            // I don't like this, but at least it's somewhat contained:
            const cds::Atom* atomCA  = this->FindAtom("CA");
            const cds::Atom* atomC   = this->FindAtom("C");
            const cds::Atom* atomO   = this->FindAtom("O");
            cds::Coordinate oxtCoord = cds::calculateCoordinateFromInternalCoords(
                atomCA->coordinate(), atomC->coordinate(), atomO->coordinate(), 120.0, 180.0, 1.25);
            this->addPdbAtom(data, residueId, "OXT", oxtCoord);
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Created new atom named OXT after " + atomO->getName() + "_" +
                          std::to_string(atomO->getNumber()));
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}
