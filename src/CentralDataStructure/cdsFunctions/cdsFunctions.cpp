#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/references.hpp"

#include <array>
#include <vector>

std::vector<int> cds::serializedNumberVector(size_t count)
{
    std::vector<int> result;
    result.reserve(count);
    for (size_t n = 0; n < count; n++)
    {
        result.push_back(n + 1);
    }
    return result;
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

cds::Sphere cds::coordinateWithRadius(Atom* atom)
{
    auto element = atom->cachedElement();
    return {MolecularMetadata::vanDerWaalsRadius(element), atom->coordinate()};
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

std::vector<std::string> cds::atomElements(const std::vector<Atom*>& atoms)
{
    std::vector<std::string> result;
    result.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        result.push_back(atom->getElement());
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

std::vector<double> cds::atomRadii(const std::vector<cds::Atom*>& atoms)
{
    std::vector<double> radii;
    radii.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        auto element = atom->cachedElement();
        radii.push_back(MolecularMetadata::vanDerWaalsRadius(element));
    }
    return radii;
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

std::vector<cds::Sphere> cds::atomCoordinatesWithRadii(const std::vector<Atom*>& atoms)
{
    std::vector<Sphere> spheres;
    spheres.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        spheres.push_back(coordinateWithRadius(atom));
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
        result.push_back(residue->getStringId());
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

std::vector<bool> cds::residuesHaveAllExpectedAtoms(const std::vector<Residue*>& residues)
{
    std::vector<bool> result;
    result.reserve(residues.size());
    for (auto& residue : residues)
    {
        result.push_back(residue->hasAllExpectedAtoms());
    }
    return result;
}

std::array<cds::Atom*, 2> cds::bondedAtomPair(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB)
{
    for (auto& a : atomsA)
    {
        for (auto& b : atomsB)
        {
            if (a->isNeighbor(b))
            {
                return {a, b};
            }
        }
    }
    return {nullptr, nullptr};
}

std::vector<bool> cds::atomsBondedTo(Atom* origin, std::vector<Atom*>& atoms)
{
    std::vector<bool> result;
    for (size_t n = 0; n < atoms.size(); n++)
    {
        result.push_back((atoms[n] == origin) || origin->isNeighbor(atoms[n]));
    }
    return result;
}
