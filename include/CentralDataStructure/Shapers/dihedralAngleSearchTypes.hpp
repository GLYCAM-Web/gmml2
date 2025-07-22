#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCHTYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLESEARCHTYPES_HPP

#include "include/CentralDataStructure/Shapers/dihedralShapeTypes.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <array>
#include <functional>
#include <vector>

namespace gmml
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

    struct DihedralIndices
    {
        std::array<size_t, 4> atoms;
        std::vector<size_t> movingAtoms;
        size_t currentMetadataIndex;
    };

    typedef std::function<std::vector<double>(const DihedralAngleData&, double, double)> SearchAngles;
    typedef std::function<double(const assembly::Bounds&)> SearchOverlap;

    struct AngleSearchSettings
    {
        double deviation;
        SearchAngles angles;
    };
} // namespace gmml
#endif
