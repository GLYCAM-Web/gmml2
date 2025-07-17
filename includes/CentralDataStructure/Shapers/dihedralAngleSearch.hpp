#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/Assembly/assemblyTypes.hpp"

#include <array>
#include <functional>
#include <stdexcept>
#include <vector>

namespace cds
{
    using GlycamMetadata::AngleDeviation;
    using GlycamMetadata::AngleLimit;
    using GlycamMetadata::AngleStd;
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::DihedralAngleDataTable;

    template<typename T>
    T onAngleDeviation(std::function<T(const AngleLimit&)>& onLimit, std::function<T(const AngleStd&)>& onStd,
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
            throw std::runtime_error("unhandled angle deviation in cds::onAngleDeviation");
        }
    };

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);
    OverlapState wiggleUsingRotamers(cds::SearchOverlap searchOverlap, SearchAngles searchAngles,
                                     const DihedralAngleDataTable& metadataTable, const assembly::Graph& graph,
                                     const assembly::Bounds& bounds, const std::vector<size_t>& movingAtoms,
                                     const DihedralCoordinates coordinates, const std::vector<size_t>& indices,
                                     const std::vector<size_t>& rotamers, const AngleSearchPreference& preference);
    assembly::Bounds simpleWiggleCurrentRotamers(
        const DihedralAngleDataTable& metadataTable, const MolecularMetadata::PotentialTable& potential,
        double overlapTolerance, SearchAngles searchAngles, std::vector<RotatableDihedral>& dihedrals,
        const std::vector<std::vector<size_t>>& metadata, const std::vector<AngleSearchPreference>& preference,
        const GraphObjects& objects, const assembly::Graph& graph, const assembly::Selection& selection,
        const assembly::Bounds& bounds, const std::vector<std::array<std::vector<bool>, 2>> residueAtomsCloseToEdge);
    std::vector<double> evenlySpacedAngles(double preference, double lowerDeviation, double upperDeviation,
                                           double increment);
    std::vector<AngleSearchPreference> angleSearchPreference(double deviation,
                                                             const ResidueLinkageShapePreference& preference);
    std::vector<std::vector<AngleSearchPreference>> angleSearchPreference(double deviation,
                                                                          const GlycanShapePreference& preferences);

    constexpr auto defaultSearchAngles = [](const DihedralAngleData& metadata, double preference, double deviation)
    {
        double increment                                              = 1.0;
        std::function<std::vector<double>(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
        {
            return cds::evenlySpacedAngles(preference, dev.lowerDeviationLimit, dev.upperDeviationLimit, increment);
        };
        std::function<std::vector<double>(const AngleStd&)> onStd = [&](const AngleStd& dev)
        {
            return cds::evenlySpacedAngles(preference, deviation * dev.lowerDeviationStd,
                                           deviation * dev.upperDeviationStd, increment);
        };
        return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
    };
    const cds::AngleSearchSettings defaultSearchSettings = {1.0, defaultSearchAngles};
} // namespace cds
#endif
