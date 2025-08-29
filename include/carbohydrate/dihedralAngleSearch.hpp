#ifndef INCLUDE_CARBOHYDRATE_DIHEDRALANGLESEARCH_HPP
#define INCLUDE_CARBOHYDRATE_DIHEDRALANGLESEARCH_HPP

#include "include/assembly/assemblyOverlap.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/metadata/elements.hpp"

#include <array>
#include <functional>
#include <stdexcept>
#include <vector>

namespace gmml
{
    template<typename T>
    T onAngleDeviation(
        std::function<T(const AngleLimit&)>& onLimit,
        std::function<T(const AngleStd&)>& onStd,
        const AngleDeviation& deviation)
    {
        if (std::holds_alternative<AngleLimit>(deviation))
        {
            return onLimit(std::get<AngleLimit>(deviation));
        }
        else if (std::holds_alternative<AngleStd>(deviation))
        {
            return onStd(std::get<AngleStd>(deviation));
        }
        else
        {
            throw std::runtime_error("unhandled angle deviation in onAngleDeviation");
        }
    };

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);

    OverlapState wiggleUsingRotamers(
        SearchOverlap searchOverlap,
        SearchAngles searchAngles,
        uint halfIntervalSearches,
        const DihedralAngleDataTable& metadataTable,
        const assembly::Graph& graph,
        const assembly::Bounds& bounds,
        const std::vector<size_t>& movingAtoms,
        const DihedralCoordinates coordinates,
        const std::vector<size_t>& indices,
        const std::vector<size_t>& rotamers,
        const AngleSearchPreference& preference);

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
        const std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge);

    std::vector<AngleSearchPreference> angleSearchPreference(
        double deviation, const ResidueLinkageShapePreference& preference);

    std::vector<std::vector<AngleSearchPreference>> angleSearchPreference(
        double deviation, const GlycanShapePreference& preferences);

    constexpr auto defaultSearchAngles =
        [](const DihedralAngleData& metadata, double preference, double deviation, uint halfIntervalSearches)
    {
        double increment = 1.0;
        std::function<AngleSpacing(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
        {
            return AngleSpacing {
                preference, dev.lowerDeviationLimit, dev.upperDeviationLimit, increment, halfIntervalSearches};
        };
        std::function<AngleSpacing(const AngleStd&)> onStd = [&](const AngleStd& dev)
        {
            return AngleSpacing {
                preference,
                deviation * dev.lowerDeviationStd,
                deviation * dev.upperDeviationStd,
                increment,
                halfIntervalSearches};
        };
        return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
    };
    const AngleSearchSettings defaultSearchSettings = {3.0, 2, defaultSearchAngles};
} // namespace gmml
#endif
