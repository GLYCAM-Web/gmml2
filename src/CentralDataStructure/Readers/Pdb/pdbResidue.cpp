#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace
#include <sstream>
#include <string>

using pdb::PdbResidue;

cds::Atom* PdbResidue::addPdbAtom(const std::string& line)
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
    atomData.recordNames.push_back(codeUtils::RemoveWhiteSpace(line.substr(0, 6)));
    atomData.occupancies.push_back(occupancy);
    atomData.temperatureFactors.push_back(temperatureFactor);
    cds::Atom* atom = this->addAtom(std::make_unique<cds::Atom>());
    readAtom(atom, line);
    return atom;
}

cds::Atom* PdbResidue::addPdbAtom(const std::string& name, const cds::Coordinate& c)
{
    atomData.recordNames.push_back("ATOM");
    atomData.occupancies.push_back(1.0);
    atomData.temperatureFactors.push_back(0.0);
    return this->addAtom(std::make_unique<cds::Atom>(name, c));
}

void PdbResidue::deletePdbAtom(cds::Atom* atom)
{
    if (atom != nullptr)
    {
        const std::vector<cds::Atom*> atoms = this->getAtoms();
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Deleting atom with id: " + atom->getName() + "_" + std::to_string(atom->getNumber()));
        size_t index = codeUtils::indexOf(atoms, atom);
        atomData.recordNames.erase(atomData.recordNames.begin() + index);
        atomData.occupancies.erase(atomData.occupancies.begin() + index);
        atomData.temperatureFactors.erase(atomData.temperatureFactors.begin() + index);
        this->deleteAtom(atom);
    }
}

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidue::PdbResidue(std::stringstream& singleResidueSecion, std::string firstLine)
{
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
                this->addPdbAtom(line);
            }
            else if (firstFoundAlternativeLocationIndicator.empty())
            { // Sometimes first atom has one location, but later atoms have alternatives. Just set and use the first
              // one (normally "A")
                firstFoundAlternativeLocationIndicator = id.getAlternativeLocation();
                this->addPdbAtom(line);
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Skipping line with alternative location: " + line);
            }
        }
    }
    this->SetType(this->determineType(this->getName()));
}

PdbResidue::PdbResidue(const std::string residueName, const PdbResidue* referenceResidue)
    : cds::Residue(residueName, referenceResidue)
{ // should instead call copy constructor and then rename with residueName?
    this->setInsertionCode(referenceResidue->getInsertionCode());
    this->setChainId(referenceResidue->getChainId());
    this->SetType(this->determineType(this->getName()));
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
pdb::ResidueId PdbResidue::getId() const
{
    ResidueId temp(this->getName(), std::to_string(this->getNumber()), this->getInsertionCode(), this->getChainId());
    return temp;
}

const std::string PdbResidue::getNumberAndInsertionCode() const
{
    return std::to_string(this->getNumber()) + this->getInsertionCode();
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void PdbResidue::modifyNTerminal(const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + this->printId());
    if (type == "NH3+")
    {
        this->deletePdbAtom(this->FindAtom("H"));
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void PdbResidue::modifyCTerminal(const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + this->printId());
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
            this->addPdbAtom("OXT", oxtCoord);
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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidue::Print(std::ostream& out) const
{
    out << "pdb::Residue : " << this->printId() << std::endl;
    for (auto& atom : this->getAtoms())
    {
        auto coord = atom->coordinate();
        out << "    atom : " << atom->getName() << "_" << atom->getNumber() << " X: " << coord.GetX()
            << " Y: " << coord.GetY() << " Z: " << coord.GetZ() << "\n";
    }
}
