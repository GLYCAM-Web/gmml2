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
    void setIntersectingIndices(std::vector<size_t>& result, cds::Sphere sphere, const std::vector<cds::Sphere>& coords,
                                const std::vector<size_t>& indices)
    {
        for (size_t index : indices)
        {
            auto& a = coords[index];
            if (cds::spheresOverlap(constants::overlapTolerance, sphere, a))
            {
                result.push_back(index);
            }
        }
    }

    void setNonIgnoredIndices(std::vector<size_t>& result, const std::vector<size_t>& indices,
                              const std::vector<bool>& ignored)
    {
        if (ignored.size() != indices.size())
        {
            throw std::runtime_error("panic");
        }
        for (size_t n = 0; n < indices.size(); n++)
        {
            if (!ignored[n])
            {
                result.push_back(indices[n]);
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
        return {coordinates,  boundingSpheres, codeUtils::indexVector(residueAtoms),
                residueAtoms, input.weights,   firstResidueBondedAtoms};
    }
} // namespace

cds::Overlap cds::CountOverlappingAtoms(const ResidueAtomOverlapInputReference& input,
                                        const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB)
{
    std::vector<size_t> indicesA;
    std::vector<size_t> indicesB;
    double tolerance = constants::overlapTolerance;
    auto properties  = OverlapProperties {constants::clashWeightBase, tolerance};
    Overlap overlap {0.0, 0.0};
    for (size_t n = 0; n < residuesA.size(); n++)
    {
        size_t aIndex = residuesA[n];
        auto& sphereA = input.boundingSpheres[aIndex];
        for (size_t k = 0; k < residuesB.size(); k++)
        {
            size_t bIndex = residuesB[k];
            auto& sphereB = input.boundingSpheres[bIndex];
            auto& atomsA  = input.residueAtoms[aIndex];
            auto& atomsB  = input.residueAtoms[bIndex];
            double weight = input.residueWeights[aIndex] * input.residueWeights[bIndex];
            indicesA.clear();
            indicesB.clear();
            if ((n == 0) && (k == 0))
            {
                setNonIgnoredIndices(indicesA, atomsA, input.firstResidueBondedAtoms[0]);
                setNonIgnoredIndices(indicesB, atomsB, input.firstResidueBondedAtoms[1]);
            }
            else if (cds::spheresOverlap(tolerance, sphereA, sphereB))
            {
                setIntersectingIndices(indicesA, sphereB, input.atomCoordinates, atomsA);
                setIntersectingIndices(indicesB, sphereA, input.atomCoordinates, atomsB);
            }
            for (size_t an : indicesA)
            {
                for (size_t bn : indicesB)
                {
                    overlap +=
                        (overlapAmount(properties, input.atomCoordinates[an], input.atomCoordinates[bn]) * weight);
                }
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

        size_t atomOffset    = inputA.atomCoordinates.size();
        size_t residueOffset = inputA.boundingSpheres.size();
        auto residueAtomsB   = inputB.residueAtoms;
        for (size_t n = 0; n < residueAtomsB.size(); n++)
        {
            residueAtomsB[n] = codeUtils::offsetIndices(atomOffset, residueAtomsB[n]);
        }

        std::vector<Sphere> atomBounds     = codeUtils::vectorAppend(inputA.atomCoordinates, inputB.atomCoordinates);
        std::vector<Sphere> residueBounds  = codeUtils::vectorAppend(inputA.boundingSpheres, inputB.boundingSpheres);
        std::vector<double> residueWeights = codeUtils::vectorAppend(inputA.residueWeights, inputB.residueWeights);
        std::vector<std::vector<size_t>> residueAtoms = codeUtils::vectorAppend(inputA.residueAtoms, residueAtomsB);
        std::array<std::vector<bool>, 2> firstBonds = {inputA.firstResidueBondedAtoms, inputB.firstResidueBondedAtoms};

        return CountOverlappingAtoms({atomBounds, residueBounds, residueAtoms, residueWeights, firstBonds},
                                     inputA.residueIndices,
                                     codeUtils::indexVectorWithOffset(residueOffset, inputB.residueIndices));
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
