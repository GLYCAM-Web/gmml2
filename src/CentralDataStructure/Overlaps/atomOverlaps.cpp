#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"

#include <numeric>

namespace
{
    void setIntersectingCoordinates(std::vector<cds::Sphere>& result, cds::Sphere sphere,
                                    const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices)
    {
        result.clear();
        for (size_t index : indices)
        {
            auto& a = coords[index];
            if (cds::spheresOverlap(constants::overlapTolerance, sphere, a))
            {
                result.push_back(a);
            }
        }
    }

    void setNonIgnoredCoordinates(std::vector<cds::Sphere>& result, const std::vector<cds::Sphere>& coords,
                                  const std::vector<size_t>& indices, const std::vector<bool>& ignored)
    {
        if (ignored.size() != indices.size())
        {
            throw std::runtime_error("panic");
        }
        result.clear();
        for (size_t n = 0; n < indices.size(); n++)
        {
            if (!ignored[n])
            {
                result.push_back(coords[indices[n]]);
            }
        }
    }

    cds::ResidueAtomOverlapInput toOverlapInput(const cds::ResiduesWithOverlapWeight& input,
                                                const std::vector<bool>& firstResidueBondedAtoms)
    {
        auto& residues   = input.residues;
        size_t atomCount = 0;
        for (auto res : residues)
        {
            atomCount += res->atomCount();
        }
        std::vector<cds::Sphere> coordinates;
        coordinates.reserve(atomCount);
        std::vector<std::vector<size_t>> residueAtoms;
        residueAtoms.reserve(residues.size());
        size_t currentAtom = 0;
        for (auto& res : residues)
        {
            auto& atomsRef = res->getAtomsReference();
            residueAtoms.push_back(codeUtils::indexVectorWithOffset(currentAtom, atomsRef));
            for (const auto& atomPtr : atomsRef)
            {
                coordinates.push_back(cds::coordinateWithRadius(atomPtr.get()));
            }
            currentAtom += atomsRef.size();
        }
        std::vector<cds::Sphere> boundingSpheres;
        boundingSpheres.reserve(residues.size());
        std::vector<cds::Sphere> residuePoints;
        for (size_t n = 0; n < residueAtoms.size(); n++)
        {
            std::vector<size_t> indices = residueAtoms[n];
            residuePoints.clear();
            for (size_t index : indices)
            {
                residuePoints.push_back(coordinates[index]);
            }
            boundingSpheres.push_back(cds::boundingSphere(residuePoints));
        }
        return {coordinates, boundingSpheres, residueAtoms, input.weights, firstResidueBondedAtoms};
    }
} // namespace

cds::Overlap cds::CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                        const ResidueAtomOverlapInputReference& moving)
{
    std::vector<Sphere> coordsA;
    std::vector<Sphere> coordsB;
    double tolerance = constants::overlapTolerance;
    auto properties  = OverlapProperties {constants::clashWeightBase, tolerance};
    Overlap overlap {0.0, 0.0};
    for (size_t n = 0; n < mostlyFixed.boundingSpheres.size(); n++)
    {
        auto& sphereA = mostlyFixed.boundingSpheres[n];
        for (size_t k = 0; k < moving.boundingSpheres.size(); k++)
        {
            double weight = mostlyFixed.residueWeights[n] * moving.residueWeights[k];
            auto& sphereB = moving.boundingSpheres[k];
            if ((n == 0) && (k == 0))
            {
                setNonIgnoredCoordinates(coordsA, mostlyFixed.atomCoordinates, mostlyFixed.residueAtoms[0],
                                         mostlyFixed.firstResidueBondedAtoms);
                setNonIgnoredCoordinates(coordsB, moving.atomCoordinates, moving.residueAtoms[0],
                                         moving.firstResidueBondedAtoms);
                overlap += (overlapAmount(properties, coordsA, coordsB) * weight);
            }
            else if (cds::spheresOverlap(tolerance, sphereA, sphereB))
            {
                setIntersectingCoordinates(coordsA, sphereB, mostlyFixed.atomCoordinates, mostlyFixed.residueAtoms[n]);
                setIntersectingCoordinates(coordsB, sphereA, moving.atomCoordinates, moving.residueAtoms[k]);
                overlap += (overlapAmount(properties, coordsA, coordsB) * weight);
            }
        }
    }
    return overlap;
}

cds::Overlap cds::CountOverlappingAtoms(const ResiduesWithOverlapWeight& residuesA,
                                        const ResiduesWithOverlapWeight& residuesB)
{
    if (residuesA.residues.empty() || residuesB.residues.empty())
    {
        return Overlap {0.0, 0.0};
    }
    else
    {
        auto bondedAtoms = [](Atom* origin, std::vector<Atom*>& atoms)
        {
            // residues aren't bonded in certain cases, use a default if so
            return (origin == nullptr) ? std::vector<bool>(atoms.size(), false) : atomsBondedTo(origin, atoms);
        };
        auto atomsA = residuesA.residues[0]->getAtoms();
        auto atomsB = residuesB.residues[0]->getAtoms();
        auto bond   = bondedAtomPair(atomsA, atomsB);
        auto inputA = toOverlapInput(residuesA, bondedAtoms(bond[0], atomsA));
        auto inputB = toOverlapInput(residuesB, bondedAtoms(bond[1], atomsB));

        return CountOverlappingAtoms({inputA.atomCoordinates, inputA.boundingSpheres, inputA.residueAtoms,
                                      inputA.residueWeights, inputA.firstResidueBondedAtoms},
                                     {inputB.atomCoordinates, inputB.boundingSpheres, inputB.residueAtoms,
                                      inputB.residueWeights, inputB.firstResidueBondedAtoms});
    }
}

cds::Overlap cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    auto coordsA = atomCoordinatesWithRadii(atomsA);
    auto coordsB = atomCoordinatesWithRadii(atomsB);

    auto properties = OverlapProperties {constants::clashWeightBase, constants::overlapTolerance};
    Overlap overlap {0, 0.0};
    for (auto& coordA : coordsA)
    {
        for (auto& coordB : coordsB)
        {
            overlap += overlapAmount(properties, coordA, coordB);
        }
    }
    return overlap;
}
