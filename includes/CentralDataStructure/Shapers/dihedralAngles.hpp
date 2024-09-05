#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <string>

namespace cds
{
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::DihedralAngleDataVector;

    struct Bounds
    {
        double lower;
        double upper;
    };

    struct AngleWithMetadata
    {
        double value;
        double preference;
        size_t metadataIndex;
    };

    std::vector<double> evenlySpaced(double lower, double upper, double approximateIncrement);
    std::vector<double> evenlySpacedAngles(double deviation, double increment, const DihedralAngleData& metadata);
    DihedralAngleDataVector likelyMetadata(const DihedralAngleDataVector& entries);
} // namespace cds
#endif
