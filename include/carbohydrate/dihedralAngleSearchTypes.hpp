#ifndef INCLUDE_CARBOHYDRATE_DIHEDRALANGLESEARCHTYPES_HPP
#define INCLUDE_CARBOHYDRATE_DIHEDRALANGLESEARCHTYPES_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralShapeTypes.hpp"
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

    struct AngleSpacing
    {
        double preference;
        double lowerDeviation;
        double upperDeviation;
        double initialIncrement;
        uint halfIntervalSearches;
    };

    typedef std::function<AngleSpacing(const DihedralAngleData&, double, double, unsigned int)> SearchAngles;
    typedef std::function<double(const assembly::Bounds&)> SearchOverlap;

    struct AngleSearchSettings
    {
        double deviation;
        uint halfIntervalSearches;
        SearchAngles angles;
    };
} // namespace gmml
#endif
