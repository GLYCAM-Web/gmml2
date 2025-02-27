#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cmath>
#include <numeric>
#include <vector>

using GlycamMetadata::DihedralAngleData;

namespace
{
    std::vector<double> evenlySpaced(double lower, double upper, double approximateIncrement)
    {
        double range = upper - lower;
        int steps    = std::ceil(std::abs(range) / approximateIncrement);
        if (steps == 0)
        {
            // range == 0, hence lower == upper
            return {};
        }
        else
        {
            double increment = range / steps;
            std::vector<double> result;
            result.reserve(steps);
            for (int k = 0; k < steps; k++)
            {
                result.push_back(lower + k * increment);
            }
            return result;
        }
    }

    std::array<std::vector<size_t>, 2> branchedResidueSets(const std::vector<bool>& residueMoving,
                                                           const std::vector<size_t>& setA,
                                                           const std::vector<size_t>& setB)
    {
        std::vector<size_t> newSetA = setB;
        std::vector<size_t> newSetB = {setA[0]};
        for (size_t n = 1; n < setA.size(); n++)
        {
            size_t index = setA[n];
            if (residueMoving[index])
            {
                newSetB.push_back(index);
            }
            else
            {
                newSetA.push_back(index);
            }
        }
        return {newSetA, newSetB};
    }

    void applyMatrix(const assembly::Graph& graph, const assembly::Bounds& initial,
                     const std::vector<size_t>& movingAtoms, const cds::RotationMatrix& matrix,
                     assembly::Bounds& bounds)
    {
        for (size_t n : movingAtoms)
        {
            bounds.atoms[n].center = matrix * initial.atoms[n].center;
        }
        assembly::updateBoundsContainingAtoms(graph, bounds, movingAtoms);
    };

    cds::AngleOverlap WiggleAnglesOverlaps(cds::SearchOverlap searchOverlap, const assembly::Graph& graph,
                                           const assembly::Bounds& initialBounds,
                                           const std::vector<size_t>& movingAtoms,
                                           const cds::DihedralCoordinates dihedral, size_t metadataIndex,
                                           double anglePreference, std::vector<double> angles)
    {
        assembly::Bounds bounds = initialBounds;
        std::vector<cds::AngleOverlap> results;
        for (double angle : angles)
        {
            cds::RotationMatrix matrix = rotationTo(dihedral, constants::toRadians(angle));
            applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
            cds::Overlap overlaps = searchOverlap(bounds);

            cds::AngleOverlap current {
                overlaps, cds::AngleWithMetadata {angle, anglePreference, metadataIndex}
            };
            if (overlaps.count <= 0.0)
            {
                // requires that the angles are sorted in order of closest to preference
                return current;
            }
            else
            {
                results.push_back(current);
            }
        }
        size_t index = bestOverlapResultIndex(results);
        return results[index];
    }
} // namespace

size_t cds::bestOverlapResultIndex(const std::vector<AngleOverlap>& results)
{
    auto differenceFromPreference = [](AngleOverlap& a)
    {
        return std::abs(a.angle.value - a.angle.preference);
    };
    size_t bestIndex = 0;
    for (size_t n = 1; n < results.size(); n++)
    {
        auto a            = results[n];
        auto b            = results[bestIndex];
        int comp          = compareOverlaps(a.overlaps, b.overlaps);
        bool sameMetadata = (a.angle.metadataIndex == b.angle.metadataIndex);
        if ((comp < 0) || (sameMetadata && (comp == 0) && differenceFromPreference(a) < differenceFromPreference(b)))
        {
            bestIndex = n;
        }
    }
    return bestIndex;
}

cds::OverlapState cds::wiggleUsingRotamers(cds::SearchOverlap searchOverlap, SearchAngles searchAngles,
                                           const assembly::Graph& graph, const assembly::Bounds& initialBounds,
                                           const std::vector<size_t>& movingAtoms,
                                           const cds::DihedralCoordinates coordinates,
                                           const std::vector<size_t>& indices, const DihedralAngleDataVector& rotamers,
                                           const AngleSearchPreference& preference)
{
    auto resultState = [&](const AngleOverlap best)
    {
        assembly::Bounds bounds    = initialBounds;
        cds::RotationMatrix matrix = rotationTo(coordinates, constants::toRadians(best.angle.value));
        applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
        return OverlapState {best.overlaps, best.angle, bounds};
    };
    std::vector<AngleOverlap> results;
    for (size_t n : preference.metadataOrder)
    {
        double angle     = preference.angles[n];
        double deviation = preference.deviation;
        AngleOverlap best =
            WiggleAnglesOverlaps(searchOverlap, graph, initialBounds, movingAtoms, coordinates, indices[n],
                                 preference.angles[n], searchAngles(rotamers[n], angle, deviation));
        // found something with no overlaps
        // if metadata and angles are sorted in order of preference, we can quit here
        if (best.overlaps.count <= 0.0)
        {
            return resultState(best);
        }
        results.push_back(best);
    }
    size_t index = bestOverlapResultIndex(results);
    return resultState(results[index]);
}

void cds::simpleWiggleCurrentRotamers(const MolecularMetadata::PotentialTable& potential, double overlapTolerance,
                                      SearchAngles searchAngles, std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata,
                                      const std::vector<AngleSearchPreference>& preference, const assembly::Graph graph,
                                      const GraphIndexData& indices,
                                      const std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge,
                                      const std::array<std::vector<Residue*>, 2>& residueSets)
{
    assembly::Bounds bounds = toAssemblyBounds(indices, graph);
    std::vector<bool> atomIncluded(graph.atomCount, true);
    std::vector<MolecularMetadata::Element> atomElements = cds::atomElements(indices.atoms);
    const std::vector<cds::Residue*>& residues           = indices.residues;
    std::vector<size_t> residueSetA                      = codeUtils::indicesOf(residues, residueSets[0]);
    std::vector<size_t> residueSetB                      = codeUtils::indicesOf(residues, residueSets[1]);
    std::vector<double> residueWeights(graph.residueCount, 1.0);
    auto dihedralCoords = [&](const cds::RotatableDihedral& dihedral)
    {
        auto coord = [&](size_t n)
        {
            return bounds.atoms[codeUtils::indexOf(indices.atoms, dihedral.atoms[n])].center;
        };
        return std::array<Coordinate, 4> {coord(3), coord(2), coord(1), coord(0)};
    };

    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        cds::RotatableDihedral& dihedral      = dihedrals[n];
        std::vector<size_t> movingAtoms       = codeUtils::indicesOf(indices.atoms, dihedral.movingAtoms);
        std::vector<Sphere> movingAtomSpheres = codeUtils::indicesToValues(bounds.atoms, movingAtoms);
        Sphere movingAtomBounds               = boundingSphere(movingAtomSpheres);
        std::array<Coordinate, 4> coordinates = dihedralCoords(dihedral);
        Coordinate pointA                     = coordinates[2];
        Coordinate pointB                     = coordinates[1];
        Sphere movementBounds                 = boundingSphereCenteredOnLine(movingAtomBounds, pointA, pointB);

        std::vector<bool> residueMoving(graph.residueCount, false);
        for (size_t n : movingAtoms)
        {
            residueMoving[graph.atomResidue[n]] = true;
        }
        std::array<std::vector<size_t>, 2> residueIndices =
            dihedral.isBranchingLinkage ? branchedResidueSets(residueMoving, residueSetA, residueSetB)
                                        : std::array<std::vector<size_t>, 2> {residueSetA, residueSetB};
        std::array<std::vector<size_t>, 2> residuesWithinBounds = {
            residueIndices[0],
            cds::intersectingIndices(overlapTolerance, movementBounds, bounds.residues, residueIndices[1])};

        std::vector<size_t> index {dihedral.currentMetadataIndex};
        auto searchOverlap = [&](const assembly::Bounds& bounds)
        {
            return cds::overlapVectorSum(cds::CountOverlappingAtoms(
                potential, overlapTolerance, graph, {bounds.atoms, atomElements, atomIncluded},
                {bounds.residues, residueWeights}, residueAtomsCloseToEdge, residuesWithinBounds[0],
                residuesWithinBounds[1]));
        };
        OverlapState best = wiggleUsingRotamers(searchOverlap, searchAngles, graph, bounds, movingAtoms, coordinates,
                                                index, metadata[n], preference[n]);
        bounds            = best.bounds;
    }
    for (size_t n = 0; n < graph.atomCount; n++)
    {
        indices.atoms[n]->setCoordinate(bounds.atoms[n].center);
    }
}

std::vector<double> cds::evenlySpacedAngles(double preference, double lowerDeviation, double upperDeviation,
                                            double increment)
{
    auto closerToPreference = [&preference](double a, double b)
    {
        return std::abs(a - preference) < std::abs(b - preference);
    };
    auto lowerRange            = evenlySpaced(preference - lowerDeviation, preference, increment);
    auto upperRange            = evenlySpaced(preference + upperDeviation, preference, increment);
    std::vector<double> result = codeUtils::vectorAppend(lowerRange, upperRange);
    result.push_back(preference);
    // sorted angles enable early return from angle search when 0 overlaps found
    std::sort(result.begin(), result.end(), closerToPreference);
    return result;
}

std::vector<cds::AngleSearchPreference> cds::angleSearchPreference(double deviation,
                                                                   const ResidueLinkageShapePreference& preference)
{
    std::function<std::vector<cds::AngleSearchPreference>(const ConformerShapePreference&)> onConformer =
        [&](const ConformerShapePreference& pref)
    {
        if (pref.metadataOrder.size() != 1)
        {
            throw std::runtime_error(
                "cannot construct search preference for conformer with multiple available rotamers");
        }
        std::vector<AngleSearchPreference> result;
        result.reserve(pref.angles.size());
        for (auto& a : pref.angles)
        {
            result.push_back({deviation, a, pref.metadataOrder});
        }
        return result;
    };
    std::function<std::vector<cds::AngleSearchPreference>(const PermutationShapePreference&)> onPermutation =
        [&](const PermutationShapePreference& pref)
    {
        std::vector<AngleSearchPreference> result;
        result.reserve(pref.angles.size());
        for (size_t n = 0; n < pref.angles.size(); n++)
        {
            result.push_back({deviation, pref.angles[n], pref.metadataOrder[n]});
        }
        return result;
    };
    return onResidueLinkageShapePreference(onConformer, onPermutation, preference);
}

std::vector<std::vector<cds::AngleSearchPreference>>
cds::angleSearchPreference(double deviation, const GlycanShapePreference& preferences)
{
    std::vector<std::vector<AngleSearchPreference>> result;
    result.reserve(preferences.size());
    for (auto& a : preferences)
    {
        result.push_back(angleSearchPreference(deviation, a));
    }
    return result;
}
