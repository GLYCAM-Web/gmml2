#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_DIHEDRALANGLES_HPP

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <string>
#include <random>

// Seed with a real random value, if available
static pcg_extras::seed_seq_from<std::random_device> seed_source;
// Make a random number engine
static pcg32 rng(seed_source);

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
        const DihedralAngleData metadata;
    };

    Bounds angleBounds(const DihedralAngleData& metadata);
    DihedralAngleDataVector likelyMetadata(const DihedralAngleDataVector& entries);
    std::string likelyName(const DihedralAngleDataVector& entries);
    AngleWithMetadata defaultAngle(const DihedralAngleData& entry);
    AngleWithMetadata randomDihedralAngleWithinMetadataRange(const DihedralAngleData& entry);
    AngleWithMetadata randomAngleEntryUsingMetadata(const DihedralAngleDataVector& entries);
} // namespace cds
#endif
