#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/references.hpp"

#include <array>
#include <vector>

std::vector<int> cds::serializedNumberVector(const std::vector<bool>& included)
{
    std::vector<int> result;
    result.reserve(included.size());
    int count = 0;
    for (bool b : included)
    {
        result.push_back(b ? count + 1 : 0);
        count += b;
    }
    return result;
}

std::vector<int> cds::serializedNumberVector(size_t count)
{
    return serializedNumberVector(std::vector<bool>(count, true));
}

size_t cds::atomVectorIndex(const std::vector<cds::Atom*>& atoms, cds::Atom* find)
{
    size_t index = codeUtils::indexOf(atoms, find);
    if (index == atoms.size())
    {
        throw std::runtime_error("atom missing from data in off writer data");
    }
    // index equals offset from start of vector
    return index;
}

std::vector<int> cds::atomNumbers(const std::vector<Atom*>& atoms)
{
    std::vector<int> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getNumber());
    }
    return result;
}

std::vector<std::string> cds::atomNames(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getName());
    }
    return result;
}

std::vector<std::string> cds::atomElementStrings(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getElement());
    }
    return result;
}

std::vector<MolecularMetadata::Element> cds::atomElements(const std::vector<Atom*>& atoms)
{
    std::vector<MolecularMetadata::Element> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->cachedElement());
    }
    return result;
}

std::vector<int> cds::atomAtomicNumbers(const std::vector<Atom*>& atoms)
{
    std::vector<int> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getAtomicNumber());
    }
    return result;
}

std::vector<cds::Coordinate> cds::atomCoordinates(const std::vector<cds::Atom*>& atoms)
{
    std::vector<Coordinate> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(atom->coordinate());
    }
    return coordinates;
}

std::vector<cds::Sphere> cds::atomCoordinatesWithRadii(const codeUtils::SparseVector<double>& elementRadii,
                                                       const std::vector<Atom*>& atoms)
{
    std::vector<Sphere> spheres;
    spheres.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        MolecularMetadata::Element element = atom->cachedElement();
        if (!elementRadii.hasValue[element])
        {
            std::string message = "No valid radius for element: " + std::to_string(element);
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        spheres.push_back({elementRadii.values[element], atom->coordinate()});
    }
    return spheres;
}

std::vector<cds::CoordinateReference> cds::atomCoordinateReferences(std::vector<cds::Atom*>& atoms)
{
    std::vector<CoordinateReference> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(atom->coordinateReference());
    }
    return coordinates;
}

std::vector<std::string> cds::atomTypes(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getType());
    }
    return result;
}

std::vector<double> cds::atomCharges(const std::vector<Atom*>& atoms)
{
    std::vector<double> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getCharge());
    }
    return result;
}

std::vector<bool> cds::atomVisibility(const std::vector<Atom*>& atoms)
{
    std::vector<bool> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->isVisible());
    }
    return result;
}

std::vector<int> cds::residueNumbers(const std::vector<Residue*>& residues)
{
    std::vector<int> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->getNumber());
    }
    return result;
}

std::vector<std::string> cds::residueNames(const std::vector<Residue*>& residues)
{
    std::vector<std::string> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->getName());
    }
    return result;
}

std::vector<std::string> cds::residueStringIds(const std::vector<Residue*>& residues)
{
    std::vector<std::string> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residueStringId(residue));
    }
    return result;
}

std::string cds::truncatedResidueName(const Residue* residue)
{
    return residue->getName().substr(0, 3);
}

std::vector<std::string> cds::truncatedResidueNames(const std::vector<Residue*>& residues)
{
    std::vector<std::string> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(truncatedResidueName(residue));
    }
    return result;
}

std::vector<cds::ResidueType> cds::residueTypes(const std::vector<Residue*>& residues)
{
    std::vector<ResidueType> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->GetType());
    }
    return result;
}

void cds::addBond(Atom* atom, Atom* otherAtom)
{
    atom->addNeighbor("atomicBond", otherAtom);
}
