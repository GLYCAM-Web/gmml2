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

    std::array<std::vector<size_t>, 2> atomsToCheck(const assembly::Graph& graph, const cds::AtomOverlapData& atomData,
                                                    const cds::ResidueOverlapData& residueData,
                                                    const std::vector<cds::BondedResidueOverlapInput>& bonds,
                                                    double overlapTolerance, size_t residueA, size_t residueB)
    {
        const cds::Sphere& residueBoundsA = residueData.bounds[residueA];
        const cds::Sphere& residueBoundsB = residueData.bounds[residueB];
        const std::vector<size_t>& atomsA = residueAtoms(graph, residueA);
        const std::vector<size_t>& atomsB = residueAtoms(graph, residueB);
        std::vector<size_t> indicesA;
        std::vector<size_t> indicesB;
        size_t bondIndex = findBondIndex(bonds, residueA, residueB);
        if (bondIndex < bonds.size())
        {
            const cds::BondedResidueOverlapInput& bond = bonds[bondIndex];
            bool order                                 = !(bond.residueIndices[0] == residueA);
            insertNonIgnored(indicesA, atomsA, atomData.included, bond.ignoredAtoms[order]);
            insertNonIgnored(indicesB, atomsB, atomData.included, bond.ignoredAtoms[!order]);
        }
        else if (cds::spheresOverlap(overlapTolerance, residueBoundsA, residueBoundsB))
        {
            insertIntersection(indicesA, overlapTolerance, residueBoundsB, atomData.bounds, atomData.included, atomsA);
            insertIntersection(indicesB, overlapTolerance, residueBoundsA, atomData.bounds, atomData.included, atomsB);
        }
        return {indicesA, indicesB};
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
                           const assembly::Graph& graph, const AtomOverlapData& atomData,
                           const ResidueOverlapData& residueData, const std::vector<BondedResidueOverlapInput>& bonds,
                           const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB)
{
    std::vector<cds::Overlap> result(graph.atomCount, {0.0, 0.0});
    for (size_t residueA : residuesA)
    {
        for (size_t residueB : residuesB)
        {
            double weight = residueData.weights[residueA] * residueData.weights[residueB];
            std::array<std::vector<size_t>, 2> toCheck =
                atomsToCheck(graph, atomData, residueData, bonds, properties.tolerance, residueA, residueB);
            for (size_t n : toCheck[0])
            {
                MolecularMetadata::Element elementA = atomData.elements[n];
                for (size_t k : toCheck[1])
                {
                    MolecularMetadata::Element elementB = atomData.elements[k];
                    double scale = MolecularMetadata::potentialWeight(potential, elementA, elementB);
                    Overlap overlap =
                        overlapAmount(properties.tolerance, scale, atomData.bounds[n], atomData.bounds[k]) * weight;
                    result[n] += overlap;
                    result[k] += overlap;
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
