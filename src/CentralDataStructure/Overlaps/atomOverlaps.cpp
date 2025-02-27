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

    std::array<std::vector<size_t>, 2>
    atomsToCheck(const assembly::Graph& graph, const cds::AtomOverlapData& atomData,
                 const cds::ResidueOverlapData& residueData,
                 const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge, double overlapTolerance,
                 size_t residueA, size_t residueB)
    {
        const cds::Sphere& residueBoundsA = residueData.bounds[residueA];
        const cds::Sphere& residueBoundsB = residueData.bounds[residueB];
        const std::vector<size_t>& atomsA = residueAtoms(graph, residueA);
        const std::vector<size_t>& atomsB = residueAtoms(graph, residueB);
        std::vector<size_t> indicesA;
        std::vector<size_t> indicesB;
        size_t adjacency                 = codeUtils::indexOf(graph.residues.nodes.nodeAdjacencies[residueA], residueB);
        const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residueA];
        if (adjacency < edges.size())
        {
            size_t edgeIndex                                     = edges[adjacency];
            const std::array<std::vector<bool>, 2>& ignoredAtoms = residueAtomsCloseToEdge[edgeIndex];
            bool order = !(graph.residues.edges.nodeAdjacencies[edgeIndex][0] == residueA);
            insertNonIgnored(indicesA, atomsA, atomData.included, ignoredAtoms[order]);
            insertNonIgnored(indicesB, atomsB, atomData.included, ignoredAtoms[!order]);
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
cds::CountOverlappingAtoms(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                           const assembly::Graph& graph, const AtomOverlapData& atomData,
                           const ResidueOverlapData& residueData,
                           const std::vector<std::array<std::vector<bool>, 2>>& residueAtomsCloseToEdge,
                           const std::vector<size_t>& residuesA, const std::vector<size_t>& residuesB)
{
    std::vector<cds::Overlap> result(graph.atomCount, {0.0, 0.0});
    for (size_t residueA : residuesA)
    {
        for (size_t residueB : residuesB)
        {
            double weight                              = residueData.weights[residueA] * residueData.weights[residueB];
            std::array<std::vector<size_t>, 2> toCheck = atomsToCheck(
                graph, atomData, residueData, residueAtomsCloseToEdge, overlapTolerance, residueA, residueB);
            for (size_t n : toCheck[0])
            {
                MolecularMetadata::Element elementA = atomData.elements[n];
                for (size_t k : toCheck[1])
                {
                    MolecularMetadata::Element elementB = atomData.elements[k];
                    double scale = MolecularMetadata::potentialWeight(potential, elementA, elementB);
                    Overlap overlap =
                        overlapAmount(overlapTolerance, scale, atomData.bounds[n], atomData.bounds[k]) * weight;
                    result[n] += overlap;
                    result[k] += overlap;
                }
            }
        }
    }
    return result;
}

cds::Overlap cds::CountOverlappingAtoms(double overlapTolerance, const std::vector<cds::Atom*>& atomsA,
                                        const std::vector<cds::Atom*>& atomsB)
{
    std::vector<Sphere> coordsA                       = atomCoordinatesWithRadii(atomsA);
    std::vector<Sphere> coordsB                       = atomCoordinatesWithRadii(atomsB);
    std::vector<MolecularMetadata::Element> elementsA = cds::atomElements(atomsA);
    std::vector<MolecularMetadata::Element> elementsB = cds::atomElements(atomsB);

    Overlap overlap {0, 0.0};
    for (size_t n = 0; n < atomsA.size(); n++)
    {
        for (size_t k = 0; k < atomsB.size(); k++)
        {
            double scale =
                MolecularMetadata::potentialWeight(MolecularMetadata::potentialTable(), elementsA[n], elementsB[k]);
            overlap += overlapAmount(overlapTolerance, scale, coordsA[n], coordsB[k]);
        }
    }
    return overlap;
}
