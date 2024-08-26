#include "includes/CentralDataStructure/cdsFunctions/atomCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <vector>

cds::Sphere cds::coordinateWithRadius(Atom* atom)
{
    auto element = atom->cachedElement();
    return {MolecularMetadata::vanDerWaalsRadius(element), *atom->getCoordinate()};
}

std::vector<cds::Coordinate*> cds::getCoordinatesFromAtoms(const std::vector<cds::Atom*>& atoms)
{
    std::vector<Coordinate*> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(atom->getCoordinate());
    }
    return coordinates;
}

std::vector<cds::Sphere> cds::getCoordinatesWithRadiiFromAtoms(const std::vector<cds::Atom*>& atoms)
{
    std::vector<cds::Sphere> coordinates;
    coordinates.reserve(atoms.size());
    for (auto& atom : atoms)
    {
        coordinates.push_back(coordinateWithRadius(atom));
    }
    return coordinates;
}
