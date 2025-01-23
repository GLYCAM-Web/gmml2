#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"

#include <numeric>

namespace
{
    void insertNonIgnored(std::vector<size_t>& result, const std::vector<size_t>& indices,
                          const std::vector<bool>& ignored)
    {
        if (ignored.size() != indices.size())
        {
            throw std::runtime_error("panic");
        }
        result.reserve(indices.size());
        for (size_t n = 0; n < indices.size(); n++)
        {
            if (!ignored[n])
            {
                result.push_back(indices[n]);
            }
        }
    }

    void insertIntersection(std::vector<size_t>& result, const cds::Sphere& sphere,
                            const std::vector<cds::Sphere>& bounds, const std::vector<size_t>& indices)
    {
        result.reserve(indices.size());
        for (size_t index : indices)
        {
            if (cds::spheresOverlap(constants::overlapTolerance, sphere, bounds[index]))
            {
                result.push_back(index);
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

    size_t findBondIndex(const std::vector<cds::BondedResidueOverlapInput>& bonds, size_t a, size_t b)
    {
        for (size_t n = 0; n < bonds.size(); n++)
        {
            const cds::BondedResidueOverlapInput& bond = bonds[n];
            size_t ax                                  = bond.residueIndices[0];
            size_t bx                                  = bond.residueIndices[1];
            if ((ax == a && bx == b) || (ax == b && bx == a))
            {
                return n;
            }
        }
        return bonds.size();
    }
} // namespace

void cds::insertIndicesOfIntersection(std::vector<size_t>& result, const Sphere& sphere,
                                      const std::vector<Sphere>& coords, const std::vector<size_t>& indices)
{
    result.reserve(indices.size());
    for (size_t index : indices)
    {
        auto& a = coords[index];
        if (cds::spheresOverlap(constants::overlapTolerance, sphere, a))
        {
            result.push_back(index);
        }
    }
}

std::vector<size_t> cds::intersectingIndices(const cds::Sphere& sphere, const std::vector<cds::Sphere>& coords,
                                             const std::vector<size_t>& indices)
{
    std::vector<size_t> result;
    insertIndicesOfIntersection(result, sphere, coords, indices);
    return result;
}

cds::Overlap cds::CountOverlappingAtoms(const std::vector<Sphere>& atomBounds, const std::vector<Sphere>& residueBounds,
                                        const std::vector<std::vector<size_t>>& residueAtoms,
                                        const std::vector<double>& residueWeights,
                                        const std::vector<BondedResidueOverlapInput>& bonds,
                                        const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB)
{
    std::vector<size_t> indicesA;
    indicesA.reserve(64);
    std::vector<size_t> indicesB;
    indicesB.reserve(64);
    double tolerance = constants::overlapTolerance;
    OverlapProperties properties {constants::clashWeightBase, tolerance};
    Overlap overlap {0.0, 0.0};
    for (size_t n = 0; n < residuesA.size(); n++)
    {
        size_t aIndex                = residuesA[n];
        const Sphere& residueBoundsA = residueBounds[aIndex];
        for (size_t k = 0; k < residuesB.size(); k++)
        {
            size_t bIndex                     = residuesB[k];
            const Sphere& residueBoundsB      = residueBounds[bIndex];
            const std::vector<size_t>& atomsA = residueAtoms[aIndex];
            const std::vector<size_t>& atomsB = residueAtoms[bIndex];
            double weight                     = residueWeights[aIndex] * residueWeights[bIndex];
            indicesA.clear();
            indicesB.clear();
            size_t bondIndex = findBondIndex(bonds, aIndex, bIndex);
            if (bondIndex < bonds.size())
            {
                const BondedResidueOverlapInput& bond = bonds[bondIndex];
                bool order                            = !(bond.residueIndices[0] == aIndex);
                insertNonIgnored(indicesA, atomsA, bond.ignoredAtoms[order]);
                insertNonIgnored(indicesB, atomsB, bond.ignoredAtoms[!order]);
            }
            else if (cds::spheresOverlap(tolerance, residueBoundsA, residueBoundsB))
            {
                insertIntersection(indicesA, residueBoundsB, atomBounds, atomsA);
                insertIntersection(indicesB, residueBoundsA, atomBounds, atomsB);
            }
            for (size_t n : indicesA)
            {
                for (size_t k : indicesB)
                {
                    overlap += (overlapAmount(properties, atomBounds[n], atomBounds[k]) * weight);
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
        std::vector<Atom*> atomsA       = residuesA.residues[0]->getAtoms();
        std::vector<Atom*> atomsB       = residuesB.residues[0]->getAtoms();
        std::array<Atom*, 2> bondedPair = bondedAtomPair(atomsA, atomsB);
        ResidueAtomOverlapInput inputA  = toOverlapInput(residuesA, bondedAtoms(bondedPair[0], atomsA));
        ResidueAtomOverlapInput inputB  = toOverlapInput(residuesB, bondedAtoms(bondedPair[1], atomsB));

        size_t atomOffset                              = inputA.atomCoordinates.size();
        size_t residueOffset                           = inputA.boundingSpheres.size();
        std::vector<std::vector<size_t>> residueAtomsB = inputB.residueAtoms;
        for (size_t n = 0; n < residueAtomsB.size(); n++)
        {
            residueAtomsB[n] = codeUtils::offsetIndices(atomOffset, residueAtomsB[n]);
        }
        std::vector<size_t> residueIndicesB = codeUtils::indexVectorWithOffset(residueOffset, inputB.residueIndices);
        std::vector<Sphere> atomBounds      = codeUtils::vectorAppend(inputA.atomCoordinates, inputB.atomCoordinates);
        std::vector<Sphere> residueBounds   = codeUtils::vectorAppend(inputA.boundingSpheres, inputB.boundingSpheres);
        std::vector<double> residueWeights  = codeUtils::vectorAppend(inputA.residueWeights, inputB.residueWeights);
        std::vector<std::vector<size_t>> residueAtoms = codeUtils::vectorAppend(inputA.residueAtoms, residueAtomsB);
        std::vector<BondedResidueOverlapInput> bonds;
        if (bondedPair[0] != nullptr)
        {
            bonds.push_back({
                {      inputA.residueIndices[0],             residueIndicesB[0]},
                {inputA.firstResidueBondedAtoms, inputB.firstResidueBondedAtoms}
            });
        }
        return CountOverlappingAtoms(atomBounds, residueBounds, residueAtoms, residueWeights, bonds,
                                     inputA.residueIndices, residueIndicesB);
    }
}

cds::Overlap cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    std::vector<Sphere> coordsA = atomCoordinatesWithRadii(atomsA);
    std::vector<Sphere> coordsB = atomCoordinatesWithRadii(atomsB);

    OverlapProperties properties {constants::clashWeightBase, constants::overlapTolerance};
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
