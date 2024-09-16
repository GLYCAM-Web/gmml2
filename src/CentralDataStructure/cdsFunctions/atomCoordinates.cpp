#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <vector>

cds::Sphere cds::coordinateWithRadius(Atom* atom)
{
    auto element = atom->cachedElement();
    return {MolecularMetadata::vanDerWaalsRadius(element), atom->coordinate()};
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

std::vector<cds::Coordinate*> cds::atomCoordinatePointers(std::vector<cds::Atom*>& atoms)
{
    std::vector<Coordinate*> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(atom->coordinatePointer());
    }
    return coordinates;
}

std::vector<cds::Sphere> cds::atomCoordinatesWithRadii(const std::vector<cds::Atom*>& atoms)
{
    std::vector<cds::Sphere> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(coordinateWithRadius(atom));
    }
    return coordinates;
}
