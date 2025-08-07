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
        const std::vector<DihedralIndices>& dihedrals,
        const std::vector<std::vector<size_t>>& metadata,
        const std::vector<AngleSearchPreference>& preference,
        const std::vector<Element>& atomElements,
        const assembly::Graph& graph,
        const assembly::Selection& selection,
        const assembly::Bounds& initialBounds,
        const std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge);

    std::vector<double> evenlySpacedAngles(
        double preference, double lowerDeviation, double upperDeviation, double increment);

    std::vector<AngleSearchPreference> angleSearchPreference(
        double deviation, const ResidueLinkageShapePreference& preference);

    std::vector<std::vector<AngleSearchPreference>> angleSearchPreference(
        double deviation, const GlycanShapePreference& preferences);

    constexpr auto defaultSearchAngles = [](const DihedralAngleData& metadata, double preference, double deviation)
    {
        double increment = 1.0;
        std::function<std::vector<double>(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
        { return evenlySpacedAngles(preference, dev.lowerDeviationLimit, dev.upperDeviationLimit, increment); };
        std::function<std::vector<double>(const AngleStd&)> onStd = [&](const AngleStd& dev)
        {
            return evenlySpacedAngles(
                preference, deviation * dev.lowerDeviationStd, deviation * dev.upperDeviationStd, increment);
        };
        return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
    };
    const AngleSearchSettings defaultSearchSettings = {1.0, defaultSearchAngles};
} // namespace gmml
#endif
