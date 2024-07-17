#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"

#include <numeric>

namespace
{

    unsigned int coordinateOverlaps(const std::vector<cds::Sphere>& coordsA, const std::vector<cds::Sphere>& coordsB)
    {
        unsigned int count = 0;
        for (auto& a : coordsA)
        {
            for (auto& b : coordsB)
            {
                count += cds::spheresOverlap(constants::overlapTolerance, a, b);
            }
        }
        return count;
    }

    void setIntersectingCoordinates(std::vector<cds::Sphere>& result, cds::Sphere sphere,
                                    const std::vector<cds::Sphere>& coords, const std::pair<size_t, size_t>& range)
    {
        result.clear();
        for (size_t n = range.first; n < range.second; n++)
        {
            auto& a = coords[n];
            if (cds::spheresOverlap(constants::overlapTolerance, sphere, a))
            {
                result.push_back(a);
            }
        }
    }
} // namespace

cds::ResidueAtomOverlapInput cds::toOverlapInput(const std::vector<Residue*>& residues)
{
    size_t atomCount = 0;
    for (auto res : residues)
    {
        atomCount += res->atomCount();
    }
    std::vector<Sphere> coordinates;
    coordinates.reserve(atomCount);
    std::vector<std::pair<size_t, size_t>> residueAtoms;
    residueAtoms.reserve(residues.size());
    size_t currentAtom = 0;
    for (auto& res : residues)
    {
        size_t startAtom = currentAtom;
        for (const auto& atomPtr : res->getAtomsReference())
        {
            coordinates.push_back(coordinateWithRadius(atomPtr.get()));
        }
        currentAtom += res->atomCount();
        residueAtoms.push_back({startAtom, currentAtom});
    }
    std::vector<Sphere> boundingSpheres;
    boundingSpheres.reserve(residues.size());
    std::vector<Sphere> residuePoints;
    for (size_t n = 0; n < residueAtoms.size(); n++)
    {
        auto range = residueAtoms[n];
        residuePoints.clear();
        residuePoints.insert(residuePoints.end(), coordinates.begin() + range.first,
                             coordinates.begin() + range.second);
        boundingSpheres.push_back(boundingSphere(residuePoints));
    }
    return {coordinates, boundingSpheres, residueAtoms};
}

unsigned int cds::CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                        const ResidueAtomOverlapInputReference& moving)
{
    std::vector<Sphere> coordsA;
    std::vector<Sphere> coordsB;
    unsigned int overlapCount = 0;
    for (size_t n = 0; n < mostlyFixed.boundingSpheres.size(); n++)
    {
        auto& sphereA = mostlyFixed.boundingSpheres[n];
        for (size_t k = 0; k < moving.boundingSpheres.size(); k++)
        {
            auto& sphereB = moving.boundingSpheres[k];
            if ((n > 0 || k > 0) && cds::spheresOverlap(constants::overlapTolerance, sphereA, sphereB))
            {
                setIntersectingCoordinates(coordsA, sphereB, mostlyFixed.atomCoordinates, mostlyFixed.residueAtoms[n]);
                setIntersectingCoordinates(coordsB, sphereA, moving.atomCoordinates, moving.residueAtoms[k]);
                overlapCount += coordinateOverlaps(coordsA, coordsB);
            }
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingAtoms(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB)
{
    auto inputA = toOverlapInput(residuesA);
    auto inputB = toOverlapInput(residuesB);

    return CountOverlappingAtoms({inputA.atomCoordinates, inputA.boundingSpheres, inputA.residueAtoms},
                                 {inputB.atomCoordinates, inputB.boundingSpheres, inputB.residueAtoms});
}

unsigned int cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    auto coordsA = getCoordinatesWithRadiiFromAtoms(atomsA);
    auto coordsB = getCoordinatesWithRadiiFromAtoms(atomsB);

    unsigned int overlapCount = 0;
    for (auto& coordA : coordsA)
    {
        for (auto& coordB : coordsB)
        {
            overlapCount += cds::spheresOverlap(constants::overlapTolerance, coordA, coordB);
        }
    }
    return overlapCount;
}
