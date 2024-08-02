#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <string>

namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

    struct Bounds
    {
        double lower;
        double upper;
    };

    struct AngleWithMetadata
    {
        double value;
        double defaultAngle;
        size_t metadataIndex;
    };

    Bounds angleBounds(const DihedralAngleData& metadata);
    DihedralAngleDataVector likelyMetadata(const DihedralAngleDataVector& entries);
} // namespace cds
#endif
