#include "include/carbohydrate/dihedralAngleSearch.hpp"

#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyOverlap.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/matrix.hpp"
#include "include/geometry/orientation.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/metadata/elements.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"

#include <cmath>
#include <numeric>
#include <vector>

namespace gmml
{
    namespace
    {
        void applyMatrix(
            const assembly::Graph& graph,
            const assembly::Bounds& initial,
            const std::vector<size_t>& movingAtoms,
            const Matrix4x4& matrix,
            assembly::Bounds& bounds)
        {
            for (size_t n : movingAtoms)
            {
                bounds.atoms[n].center = matrix * initial.atoms[n].center;
            }
            assembly::updateBoundsContainingAtoms(graph, bounds, movingAtoms);
        };

        struct AngleWithAdjacent
        {
            double angle;
            double lowerIncrement;
            double upperIncrement;
        };

        AngleOverlap findLeastOverlapAngle(
            SearchOverlap searchOverlap,
            const assembly::Graph& graph,
            const assembly::Bounds& initialBounds,
            const std::vector<size_t>& movingAtoms,
            const DihedralCoordinates dihedral,
            size_t metadataIndex,
            const AngleSpacing& spacing)
        {
            double preference = spacing.preference;
            assembly::Bounds bounds = initialBounds;
            // lambda reuses and mutates bounds variable to avoid repeated memory allocation
            // without parallel execution it's fine
            auto overlapAt = [&](double angle)
            {
                Matrix4x4 matrix = rotationTo(dihedral, constants::toRadians(angle));
                applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
                double overlaps = searchOverlap(bounds);

                return AngleOverlap {
                    overlaps, {angle, preference, metadataIndex}
                };
            };
            AngleOverlap def = overlapAt(preference);
            if (def.overlaps == 0)
            {
                return def;
            }

            double lowerRange = spacing.lowerDeviation;
            double upperRange = spacing.upperDeviation;
            uint lowerSteps = std::floor(lowerRange / spacing.initialIncrement);
            uint upperSteps = std::floor(upperRange / spacing.initialIncrement);
            double lowerIncrement = lowerRange / (lowerSteps + 1);
            double upperIncrement = upperRange / (upperSteps + 1);

            std::vector<AngleWithAdjacent> anglesToCheck;
            anglesToCheck.reserve(lowerSteps + upperSteps + 1);
            anglesToCheck.push_back({preference, lowerIncrement, upperIncrement});
            for (size_t n = 1; n <= lowerSteps; n++)
            {
                anglesToCheck.push_back({preference - n * lowerIncrement, lowerIncrement, lowerIncrement});
            }
            for (size_t n = 1; n <= upperSteps; n++)
            {
                anglesToCheck.push_back({preference + n * upperIncrement, upperIncrement, upperIncrement});
            }

            std::sort(
                anglesToCheck.begin(),
                anglesToCheck.end(),
                [&](const AngleWithAdjacent& a, const AngleWithAdjacent& b)
                { return std::abs(a.angle - preference) < (b.angle - preference); });

            std::vector<AngleOverlap> results {def};

            for (size_t n = 1; n < anglesToCheck.size(); n++)
            {
                double angle = anglesToCheck[n].angle;
                AngleOverlap current = overlapAt(angle);
                results.push_back(current);
                if (current.overlaps <= 0.0)
                {
                    break;
                }
            }
            size_t index = bestOverlapResultIndex(results);
            AngleWithAdjacent lookBy = anglesToCheck[index];
            results = {results[index]};
            for (size_t k = 0; k < spacing.halfIntervalSearches; k++)
            {
                double lowerInc = lookBy.lowerIncrement * 0.5;
                double upperInc = lookBy.upperIncrement * 0.5;
                double lower = lookBy.angle - lowerInc;
                double upper = lookBy.angle + upperInc;
                std::vector<AngleWithAdjacent> localAngles {
                    {lower, lowerInc, lowerInc},
                    {upper, upperInc, upperInc}
                };
                std::vector<AngleOverlap> local {overlapAt(lower), overlapAt(upper)};
                size_t bestLocal = bestOverlapResultIndex(local);
                results.push_back(local[bestLocal]);
                lookBy = localAngles[bestLocal];
            }
            size_t bestIndex = bestOverlapResultIndex(results);
            return results[bestIndex];
        }
    } // namespace

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results)
    {
        auto differenceFromPreference = [](AngleOverlap& a) { return std::abs(a.angle.value - a.angle.preference); };
        size_t bestIndex = 0;
        for (size_t n = 1; n < results.size(); n++)
        {
            auto a = results[n];
            auto b = results[bestIndex];
            int comp = compareOverlaps(a.overlaps, b.overlaps);
            bool sameMetadata = (a.angle.metadataIndex == b.angle.metadataIndex);
            if ((comp < 0) ||
                (sameMetadata && (comp == 0) && differenceFromPreference(a) < differenceFromPreference(b)))
            {
                bestIndex = n;
            }
        }
        return bestIndex;
    }

    OverlapState wiggleUsingRotamers(
        SearchOverlap searchOverlap,
        SearchAngles searchAngles,
        uint halfIntervalSearches,
        const DihedralAngleDataTable& metadataTable,
        const assembly::Graph& graph,
        const assembly::Bounds& initialBounds,
        const std::vector<size_t>& movingAtoms,
        const DihedralCoordinates coordinates,
        const std::vector<size_t>& indices,
        const std::vector<size_t>& rotamers,
        const AngleSearchPreference& preference)
    {
        auto resultState = [&](const AngleOverlap best)
        {
            assembly::Bounds bounds = initialBounds;
            Matrix4x4 matrix = rotationTo(coordinates, constants::toRadians(best.angle.value));
            applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
            return OverlapState {best.overlaps, best.angle, bounds};
        };
        std::vector<AngleOverlap> results;
        for (size_t n : preference.metadataOrder)
        {
            double angle = preference.angles[n];
            double deviation = preference.deviation;
            AngleOverlap best = findLeastOverlapAngle(
                searchOverlap,
                graph,
                initialBounds,
                movingAtoms,
                coordinates,
                indices[n],
                searchAngles(metadataTable.entries[rotamers[n]], angle, deviation, halfIntervalSearches));
            // found something with no overlaps
            // if metadata and angles are sorted in order of preference, we can quit here
            if (best.overlaps <= 0.0)
            {
                return resultState(best);
            }
            results.push_back(best);
        }
        size_t index = bestOverlapResultIndex(results);
        return resultState(results[index]);
    }

    assembly::Bounds simpleWiggleCurrentRotamers(
        const DihedralAngleDataTable& metadataTable,
        const PotentialTable& potential,
        SearchAngles searchAngles,
        uint halfIntervalSearches,
        const std::vector<DihedralIndices>& dihedrals,
        const std::vector<std::vector<size_t>>& metadata,
        const std::vector<AngleSearchPreference>& preference,
        const std::vector<Element>& atomElements,
        const assembly::Graph& graph,
        const assembly::Selection& selection,
        const assembly::Bounds& initialBounds,
        const std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge)
    {
        assembly::Bounds bounds = initialBounds;
        auto dihedralCoords = [&](const DihedralIndices& dihedral)
        {
            auto coord = [&](size_t n) { return bounds.atoms[dihedral.atoms[n]].center; };
            return std::array<Coordinate, 4> {coord(3), coord(2), coord(1), coord(0)};
        };

        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            const DihedralIndices& dihedral = dihedrals[n];
            std::vector<bool> atomMoving = util::indicesToBools(atomCount(graph.source), dihedral.movingAtoms);
            std::array<Coordinate, 4> coordinates = dihedralCoords(dihedral);
            assembly::Selection selectionA =
                assembly::intersection(graph, selection, assembly::selectByAtoms(graph, atomMoving));
            assembly::Selection selectionB =
                assembly::intersection(graph, selection, assembly::selectByAtoms(graph, util::vectorNot(atomMoving)));

            std::vector<size_t> index {dihedral.currentMetadataIndex};

            auto searchOverlap = [&](const assembly::Bounds& bounds)
            {
                return overlapVectorSum(overlapsBetweenSelections(
                    potential, 0.0, graph, bounds, selectionA, selectionB, atomElements, residueAtomsCloseToEdge));
            };
            OverlapState best = wiggleUsingRotamers(
                searchOverlap,
                searchAngles,
                halfIntervalSearches,
                metadataTable,
                graph,
                bounds,
                dihedral.movingAtoms,
                coordinates,
                index,
                metadata[n],
                preference[n]);
            bounds = best.bounds;
        }
        return bounds;
    }

    std::vector<AngleSearchPreference> angleSearchPreference(
        double deviation, const ResidueLinkageShapePreference& preference)
    {
        std::function<std::vector<AngleSearchPreference>(const ConformerShapePreference&)> onConformer =
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
        std::function<std::vector<AngleSearchPreference>(const PermutationShapePreference&)> onPermutation =
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

    std::vector<std::vector<AngleSearchPreference>> angleSearchPreference(
        double deviation, const GlycanShapePreference& preferences)
    {
        std::vector<std::vector<AngleSearchPreference>> result;
        result.reserve(preferences.size());
        for (auto& a : preferences)
        {
            result.push_back(angleSearchPreference(deviation, a));
        }
        return result;
    }
} // namespace gmml
