#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"

#include <numeric>

namespace
{
    cds::Overlap sphereOverlap(double tolerance, const cds::Sphere& a, const cds::Sphere& b)
    {
        double cutoff = a.radius + b.radius - tolerance;
        double sqDist = cds::squaredDistance(a.center, b.center);
        double base   = constants::clashWeightBase;
        if (sqDist < cutoff * cutoff)
        {
            return cds::Overlap {1.0, base / (base + sqDist)};
        }
        else
        {
            return cds::Overlap {0.0, 0.0};
        }
    }

    cds::Overlap coordinateOverlaps(const std::vector<cds::Sphere>& coordsA, const std::vector<cds::Sphere>& coordsB)
    {
        cds::Overlap overlap {0.0, 0.0};
        for (auto& a : coordsA)
        {
            for (auto& b : coordsB)
            {
                overlap += sphereOverlap(constants::overlapTolerance, a, b);
            }
        }
        return overlap;
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

cds::ResidueAtomOverlapInput cds::toOverlapInput(const ResiduesWithOverlapWeight& input)
{
    auto& residues   = input.residues;
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
    return {coordinates, boundingSpheres, residueAtoms, input.weights};
}

cds::Overlap cds::CountOverlappingAtoms(bool ignoreNeighboringResidues,
                                        const ResidueAtomOverlapInputReference& mostlyFixed,
                                        const ResidueAtomOverlapInputReference& moving)
{
    std::vector<Sphere> coordsA;
    std::vector<Sphere> coordsB;
    Overlap overlap {0.0, 0.0};
    for (size_t n = 0; n < mostlyFixed.boundingSpheres.size(); n++)
    {
        auto& sphereA = mostlyFixed.boundingSpheres[n];
        for (size_t k = 0; k < moving.boundingSpheres.size(); k++)
        {
            auto& sphereB = moving.boundingSpheres[k];
            if (!(ignoreNeighboringResidues && (n == 0) && (k == 0)) &&
                cds::spheresOverlap(constants::overlapTolerance, sphereA, sphereB))
            {
                setIntersectingCoordinates(coordsA, sphereB, mostlyFixed.atomCoordinates, mostlyFixed.residueAtoms[n]);
                setIntersectingCoordinates(coordsB, sphereA, moving.atomCoordinates, moving.residueAtoms[k]);
                uint weight = mostlyFixed.residueWeights[n] * moving.residueWeights[k];
                overlap     += (coordinateOverlaps(coordsA, coordsB) * weight);
            }
        }
    }
    return overlap;
}

cds::Overlap cds::CountOverlappingAtoms(bool ignoreNeighboringResidues, const ResiduesWithOverlapWeight& residuesA,
                                        const ResiduesWithOverlapWeight& residuesB)
{
    auto inputA = toOverlapInput(residuesA);
    auto inputB = toOverlapInput(residuesB);

    return CountOverlappingAtoms(
        ignoreNeighboringResidues,
        {inputA.atomCoordinates, inputA.boundingSpheres, inputA.residueAtoms, inputA.residueWeights},
        {inputB.atomCoordinates, inputB.boundingSpheres, inputB.residueAtoms, inputB.residueWeights});
}

cds::Overlap cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    auto coordsA = getCoordinatesWithRadiiFromAtoms(atomsA);
    auto coordsB = getCoordinatesWithRadiiFromAtoms(atomsB);

    Overlap overlap {0, 0.0};
    for (auto& coordA : coordsA)
    {
        for (auto& coordB : coordsB)
        {
            overlap += sphereOverlap(constants::overlapTolerance, coordA, coordB);
        }
    }
    return overlap;
}
