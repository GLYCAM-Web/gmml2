#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <iostream>
#include <iomanip>

namespace
{
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

int pdb::checkShiftFromSerialNumberOverrun(const std::string& line)
{
    int shift = 0;
    if (isdigit(line[11]) && line[20] != ' ')
    {
        shift =
            codeUtils::GetSizeOfIntInString(line.substr(12)); // The shift starts when this gets overrun, not 11. Right?
        std::stringstream ss;
        ss << "Shift of size " << shift << " detected as position 12 is a digit: " << line[11] << " and position 21 >>>"
           << line[20] << "<<< isn't blank: " << line[20] << ":\n           v        v\n"
           << line << "\n";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
    }
    return shift;
}

int pdb::checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift)
{
    return codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
}

cds::Coordinate pdb::checkShiftsAndExtractCoordinate(const std::string& line)
{
    int shift       = pdb::checkShiftFromSerialNumberOverrun(line);
    int secondShift = pdb::checkSecondShiftFromResidueNumberOverrun(line, shift);
    // Coordinates etc don't get shifted by first overrun in residue number, but do by the rest.
    if (secondShift > 1)
    { // Combine the shifts, but ignore the first shift in residue sequence number.
        shift += (secondShift - 1);
    }
    return pdb::coordinateFromStrings(codeUtils::RemoveWhiteSpace(line.substr(30 + shift, 8)),
                                      codeUtils::RemoveWhiteSpace(line.substr(38 + shift, 8)),
                                      codeUtils::RemoveWhiteSpace(line.substr(46 + shift, 8)));
}

cds::Coordinate pdb::coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z)
{
    try
    {
        return {std::stod(x), std::stod(y), std::stod(z)};
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Could not convert these strings to doubles: " + x + ", " + y + ", " + z + ", ");
        throw;
    }
}

void pdb::expandLine(std::string& line, int length)
{
    int l = line.length();
    if (l < length)
    {
        int space = length - l;
        std::stringstream ss;
        ss << line << std::setw(space) << " ";
        line = ss.str();
    }
}

size_t pdb::addAtom(PdbData& data, size_t residueId, const AtomEntry& entry)
{
    cds::Atom* atom = data.objects.residues[residueId]->addAtom(std::make_unique<cds::Atom>());
    size_t atomId   = data.indices.atomCount;
    data.objects.atoms.push_back(atom);
    data.indices.atomResidue.push_back(residueId);
    data.indices.atomAlive.push_back(true);
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

size_t pdb::addAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& coordinate)
{
    double occupancy         = 1.0;
    double temperatureFactor = 0.0;
    return addAtom(data, residueId, {"ATOM", name, 1, coordinate, occupancy, temperatureFactor});
}

void pdb::deleteAtom(PdbData& data, size_t residueId, size_t atomId)
{
    if (atomId < data.indices.atomCount)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Deleting atom with id: " + data.atoms.names[atomId] + "_" +
                      std::to_string(data.atoms.numbers[atomId]));
        data.indices.atomAlive[atomId]   = false;
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
    std::vector<size_t>& order = data.molecules.residueOrder[moleculeId];
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

void pdb::addBond(PdbData& data, size_t atom1, size_t atom2)
{
    cds::addBond(data.objects.atoms[atom1], data.objects.atoms[atom2]);
    graph::addEdge(data.atomGraph, {atom1, atom2});
}

size_t pdb::findResidueAtom(const PdbData& data, size_t residueId, const std::string& atomName)
{
    for (size_t n : assembly::residueAtoms(data.indices, residueId))
    {
        if (data.atoms.names[n] == atomName)
        {
            return n;
        }
    }
    return data.indices.atomCount;
}
