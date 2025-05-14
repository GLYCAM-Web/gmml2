#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/biology.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace

#include <string>

namespace
{
    using cds::ResidueType;

    ResidueType residueType(const std::string& name)
    {
        if (codeUtils::contains(biology::proteinResidueNames, name))
        {
            return ResidueType::Protein;
        }
        // ToDo we want to figure out solvent, aglycone etc here too?.
        return ResidueType::Undefined;
    }

    MolecularMetadata::Element atomElement(const std::string& name)
    {
        if (isalpha(name.at(0))) // if first char is in the alphabet
        {
            return MolecularMetadata::toElement(name.substr(0, 1));
        }
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Did not find an element for atom named: " + name);
        return MolecularMetadata::Unknown;
    }
} // namespace

pdb::ResidueId pdb::pdbResidueId(const PdbData& data, size_t residueId)
{
    return ResidueId(data.residues.names[residueId], std::to_string(data.residues.numbers[residueId]),
                     data.residues.insertionCodes[residueId], data.residues.chainIds[residueId]);
}

std::string pdb::residueStringId(const PdbData& data, size_t residueId)
{
    return pdbResidueId(data, residueId).print();
}

size_t pdb::addPdbAtom(PdbData& data, size_t residueId, const AtomEntry& entry)
{
    cds::Atom* atom = data.objects.residues[residueId]->addAtom(std::make_unique<cds::Atom>());
    size_t atomId   = data.indices.atomCount;
    data.objects.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.indices.atomCount++;
    data.atoms.recordNames.push_back(entry.recordName);
    data.atoms.names.push_back(entry.name);
    data.atoms.elements.push_back(atomElement(entry.name));
    data.atoms.numbers.push_back(entry.number);
    data.atoms.coordinates.push_back(entry.coordinate);
    data.atoms.occupancies.push_back(entry.occupancy);
    data.atoms.temperatureFactors.push_back(entry.temperatureFactor);
    atom->setNumber(entry.number);
    atom->setName(entry.name);
    atom->setCoordinate(entry.coordinate);
    atom->setElement(data.atoms.elements[atomId]);
    size_t nodeId = graph::addNode(data.atomGraph);
    if (nodeId != atomId)
    {
        throw std::runtime_error("atom id mismatch");
    }
    return atomId;
}

size_t pdb::addPdbAtom(PdbData& data, size_t residueId, const std::string& line)
{
    AtomEntry entry = readAtom(line);
    return addPdbAtom(data, residueId, entry);
}

size_t pdb::addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& coordinate)
{
    double occupancy         = 1.0;
    double temperatureFactor = 0.0;
    return addPdbAtom(data, residueId, {"ATOM", name, 1, coordinate, occupancy, temperatureFactor});
}

void pdb::deletePdbAtom(PdbData& data, size_t residueId, size_t atomId)
{
    if (atomId < data.indices.atomCount)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Deleting atom with id: " + data.atoms.names[atomId] + "_" +
                      std::to_string(data.atoms.numbers[atomId]));
        data.atomGraph.nodeAlive[atomId] = false;
        cds::Atom* atom                  = data.objects.atoms[atomId];
        data.objects.residues[residueId]->deleteAtom(atom);
    }
}

size_t pdb::addResidue(PdbData& data, size_t moleculeId, size_t position, const ResidueEntry& entry)
{
    size_t residueId = data.indices.residueCount;
    cds::Residue* residue =
        data.objects.molecules[moleculeId]->insertNewResidue(std::make_unique<cds::Residue>(), position);
    std::vector<size_t>& order = data.moleculeResidueOrder[moleculeId];
    order.insert(order.begin() + position, residueId);
    data.indices.residueMolecule.push_back(moleculeId);
    data.indices.residueCount++;
    data.objects.residues.push_back(residue);
    data.residues.names.push_back(entry.name);
    data.residues.types.push_back(entry.type);
    data.residues.numbers.push_back(entry.number);
    data.residues.insertionCodes.push_back(entry.insertionCode);
    data.residues.chainIds.push_back(entry.chainId);
    data.residues.isCTerminal.push_back(false);
    data.residues.isNTerminal.push_back(false);
    data.residues.hasTerCard.push_back(entry.hasTerCard);
    residue->setName(entry.name);
    residue->SetType(entry.type);
    residue->setNumber(entry.number);
    return residueId;
}

size_t pdb::readResidue(PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine)
{
    ResidueId resId(firstLine);
    size_t position                                    = data.moleculeResidueOrder[moleculeId].size();
    std::string name                                   = resId.getName();
    size_t residueId                                   = addResidue(data, moleculeId, position,
                                                                    {name, residueType(name), uint(std::stoi(resId.getNumber())),
                                                                     resId.getInsertionCode(), resId.getChainId(), false});
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
    return residueId;
}

size_t pdb::readResidue(PdbData& data, size_t moleculeId, const std::string& name, cds::ResidueType type,
                        bool hasTerCard, size_t referenceResidue)
{
    size_t position = codeUtils::indexOf(data.moleculeResidueOrder[moleculeId], referenceResidue) + 1;
    return addResidue(data, moleculeId, position,
                      {name, type, data.residues.numbers[referenceResidue] + 1,
                       data.residues.insertionCodes[referenceResidue], data.residues.chainIds[referenceResidue],
                       hasTerCard});
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string pdb::getNumberAndInsertionCode(const PdbData& data, size_t residueId)
{
    return std::to_string(data.residues.numbers[residueId]) + data.residues.insertionCodes[residueId];
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void pdb::modifyNTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "NH3+")
    {
        size_t atomId = findResidueAtom(data, residueId, "H");
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
        return data.atoms.coordinates[n];
    };
    gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + pdbResidueId(data, residueId).print());
    if (type == "CO2-")
    {
        size_t atomOXT = findAtom("OXT");
        if (atomOXT == data.indices.atomCount)
        {
            // I don't like this, but at least it's somewhat contained:
            size_t atomCA            = findAtom("CA");
            size_t atomC             = findAtom("C");
            size_t atomO             = findAtom("O");
            cds::Coordinate oxtCoord = cds::calculateCoordinateFromInternalCoords(coord(atomCA), coord(atomC),
                                                                                  coord(atomO), 120.0, 180.0, 1.25);
            addPdbAtom(data, residueId, "OXT", oxtCoord);
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Created new atom named OXT after " + data.atoms.names[atomO] + "_" +
                          std::to_string(data.atoms.numbers[atomO]));
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}
