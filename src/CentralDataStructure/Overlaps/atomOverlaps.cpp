#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
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
                          const std::vector<bool>& includeAtom, const std::vector<bool>& localIgnore)
    {
        if (localIgnore.size() != indices.size())
        {
            throw std::runtime_error("panic");
        }
        result.reserve(indices.size());
        for (size_t n = 0; n < indices.size(); n++)
        {
            if (includeAtom[indices[n]] && !localIgnore[n])
            {
                result.push_back(indices[n]);
            }
        }
    }

    void insertIntersection(std::vector<size_t>& result, double overlapTolerance, const cds::Sphere& sphere,
                            const std::vector<cds::Sphere>& bounds, const std::vector<bool>& includeAtom,
                            const std::vector<size_t>& indices)
    {
        result.reserve(indices.size());
        for (size_t index : indices)
        {
            if (includeAtom[index] && cds::spheresOverlap(overlapTolerance, sphere, bounds[index]))
            {
                result.push_back(index);
            }
        }
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

void cds::insertIndicesOfIntersection(std::vector<size_t>& result, double overlapTolerance, const Sphere& sphere,
                                      const std::vector<Sphere>& coords, const std::vector<size_t>& indices)
{
    result.reserve(indices.size());
    for (size_t index : indices)
    {
        auto& a = coords[index];
        if (cds::spheresOverlap(overlapTolerance, sphere, a))
        {
            result.push_back(index);
        }
    }
}

std::vector<size_t> cds::intersectingIndices(double overlapTolerance, const cds::Sphere& sphere,
                                             const std::vector<cds::Sphere>& coords, const std::vector<size_t>& indices)
{
    std::vector<size_t> result;
    insertIndicesOfIntersection(result, overlapTolerance, sphere, coords, indices);
    return result;
}

std::vector<cds::Overlap>
cds::CountOverlappingAtoms(const MolecularMetadata::PotentialTable& potential, OverlapProperties properties,
                           const std::vector<Sphere>& atomBounds, const std::vector<Sphere>& residueBounds,
                           const std::vector<std::vector<size_t>>& residueAtoms,
                           const std::vector<double>& residueWeights,
                           const std::vector<MolecularMetadata::Element>& atomElements,
                           const std::vector<bool>& includedAtoms, const std::vector<BondedResidueOverlapInput>& bonds,
                           const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB)
{
    std::vector<cds::Overlap> result(atomBounds.size(), {0.0, 0.0});
    std::vector<size_t> indicesA;
    indicesA.reserve(64);
    std::vector<size_t> indicesB;
    indicesB.reserve(64);
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
                insertNonIgnored(indicesA, atomsA, includedAtoms, bond.ignoredAtoms[order]);
                insertNonIgnored(indicesB, atomsB, includedAtoms, bond.ignoredAtoms[!order]);
            }
            else if (cds::spheresOverlap(properties.tolerance, residueBoundsA, residueBoundsB))
            {
                insertIntersection(indicesA, properties.tolerance, residueBoundsB, atomBounds, includedAtoms, atomsA);
                insertIntersection(indicesB, properties.tolerance, residueBoundsA, atomBounds, includedAtoms, atomsB);
            }
            for (size_t n : indicesA)
            {
                MolecularMetadata::Element elementA = atomElements[n];
                for (size_t k : indicesB)
                {
                    MolecularMetadata::Element elementB = atomElements[k];
                    double scale    = MolecularMetadata::potentialWeight(potential, elementA, elementB);
                    Overlap overlap = overlapAmount(properties.tolerance, scale, atomBounds[n], atomBounds[k]) * weight;
                    result[n]       += overlap;
                    result[k]       += overlap;
                }
            }
        }
    }
    return result;
}

cds::Overlap cds::CountOverlappingAtoms(OverlapProperties properties, const std::vector<cds::Atom*>& atomsA,
                                        const std::vector<cds::Atom*>& atomsB)
{
    std::vector<Sphere> coordsA                       = atomCoordinatesWithRadii(atomsA);
    std::vector<Sphere> coordsB                       = atomCoordinatesWithRadii(atomsB);
    std::vector<MolecularMetadata::Element> elementsA = cds::atomElementEnums(atomsA);
    std::vector<MolecularMetadata::Element> elementsB = cds::atomElementEnums(atomsB);

    Overlap overlap {0, 0.0};
    for (size_t n = 0; n < atomsA.size(); n++)
    {
        for (size_t k = 0; k < atomsB.size(); k++)
        {
            double scale =
                MolecularMetadata::potentialWeight(MolecularMetadata::potentialTable(), elementsA[n], elementsB[k]);
            overlap += overlapAmount(properties.tolerance, scale, coordsA[n], coordsB[k]);
        }
    }
    return overlap;
}
