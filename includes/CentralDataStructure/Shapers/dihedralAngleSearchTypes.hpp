#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCHTYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCHTYPES_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/Assembly/assemblyTypes.hpp"

#include <array>
#include <vector>
#include <functional>

namespace cds
{
    struct AngleOverlap
    {
        cds::Overlap overlaps;
        AngleWithMetadata angle;
    };

    struct OverlapState
    {
        Overlap overlap;
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
    typedef std::function<cds::Overlap(const assembly::Bounds&)> SearchOverlap;

    struct AngleSearchSettings
    {
        double deviation;
        SearchAngles angles;
    };
} // namespace cds
#endif
