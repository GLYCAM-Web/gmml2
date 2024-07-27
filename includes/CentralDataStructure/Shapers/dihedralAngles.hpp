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
        DihedralAngleData metadata;
    };

    Bounds angleBounds(const DihedralAngleData& metadata);
    DihedralAngleDataVector likelyMetadata(const DihedralAngleDataVector& entries);
    std::string likelyName(const DihedralAngleDataVector& entries);
    AngleWithMetadata defaultAngle(const DihedralAngleData& entry);
} // namespace cds
#endif
