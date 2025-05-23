#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCH_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"

#include <array>
#include <vector>
#include <functional>

namespace cds
{
    struct AngleOverlap
    {
        double overlaps;
        AngleWithMetadata angle;
    };

    struct OverlapState
    {
        double overlap;
        AngleWithMetadata angle;
        assembly::Bounds bounds;
    };

    struct AngleSearchPreference
    {
        double deviation;
        std::vector<double> angles;
        std::vector<size_t> metadataOrder;
    };

    typedef std::function<std::vector<double>(const GlycamMetadata::DihedralAngleData&, double, double)> SearchAngles;
    typedef std::function<double(const assembly::Bounds&)> SearchOverlap;

    struct AngleSearchSettings
    {
        double deviation;
        SearchAngles angles;
    };

    size_t bestOverlapResultIndex(const std::vector<AngleOverlap>& results);
    OverlapState wiggleUsingRotamers(cds::SearchOverlap searchOverlap, SearchAngles searchAngles,
                                     const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                     const assembly::Graph& graph, const assembly::Bounds& bounds,
                                     const std::vector<size_t>& movingAtoms, const DihedralCoordinates coordinates,
                                     const std::vector<size_t>& indices, const std::vector<size_t>& rotamers,
                                     const AngleSearchPreference& preference);
    assembly::Bounds simpleWiggleCurrentRotamers(
        const GlycamMetadata::DihedralAngleDataTable& metadataTable, const MolecularMetadata::PotentialTable& potential,
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
} // namespace cds
#endif
