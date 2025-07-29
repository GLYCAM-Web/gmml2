#include "include/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"

#include "include/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/geometry/overlap.hpp"
#include "include/geometry/rotationMatrix.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/metadata/elements.hpp"
#include "include/structure/atomOverlaps.hpp"
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
        std::vector<double> evenlySpaced(double lower, double upper, double approximateIncrement)
        {
            double range = upper - lower;
            int steps = std::ceil(std::abs(range) / approximateIncrement);
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

        void applyMatrix(
            const assembly::Graph& graph,
            const assembly::Bounds& initial,
            const std::vector<size_t>& movingAtoms,
            const RotationMatrix& matrix,
            assembly::Bounds& bounds)
        {
            for (size_t n : movingAtoms)
            {
                bounds.atoms[n].center = matrix * initial.atoms[n].center;
            }
            assembly::updateBoundsContainingAtoms(graph, bounds, movingAtoms);
        };

        AngleOverlap WiggleAnglesOverlaps(
            SearchOverlap searchOverlap,
            const assembly::Graph& graph,
            const assembly::Bounds& initialBounds,
            const std::vector<size_t>& movingAtoms,
            const DihedralCoordinates dihedral,
            size_t metadataIndex,
            double anglePreference,
            std::vector<double> angles)
        {
            assembly::Bounds bounds = initialBounds;
            std::vector<AngleOverlap> results;
            for (double angle : angles)
            {
                RotationMatrix matrix = rotationTo(dihedral, constants::toRadians(angle));
                applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
                double overlaps = searchOverlap(bounds);

                AngleOverlap current {
                    overlaps, AngleWithMetadata {angle, anglePreference, metadataIndex}
                };
                if (overlaps <= 0.0)
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
            RotationMatrix matrix = rotationTo(coordinates, constants::toRadians(best.angle.value));
            applyMatrix(graph, initialBounds, movingAtoms, matrix, bounds);
            return OverlapState {best.overlaps, best.angle, bounds};
        };
        std::vector<AngleOverlap> results;
        for (size_t n : preference.metadataOrder)
        {
            double angle = preference.angles[n];
            double deviation = preference.deviation;
            AngleOverlap best = WiggleAnglesOverlaps(
                searchOverlap,
                graph,
                initialBounds,
                movingAtoms,
                coordinates,
                indices[n],
                preference.angles[n],
                searchAngles(metadataTable.entries[rotamers[n]], angle, deviation));
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
            std::vector<bool> atomMoving = util::indicesToBools(graph.indices.atomCount, dihedral.movingAtoms);
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

    std::vector<double> evenlySpacedAngles(
        double preference, double lowerDeviation, double upperDeviation, double increment)
    {
        auto closerToPreference = [&preference](double a, double b)
        { return std::abs(a - preference) < std::abs(b - preference); };
        auto lowerRange = evenlySpaced(preference - lowerDeviation, preference, increment);
        auto upperRange = evenlySpaced(preference + upperDeviation, preference, increment);
        std::vector<double> result = util::vectorAppend(lowerRange, upperRange);
        result.push_back(preference);
        // sorted angles enable early return from angle search when 0 overlaps found
        std::sort(result.begin(), result.end(), closerToPreference);
        return result;
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
