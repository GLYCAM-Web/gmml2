#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <string>
#include <sstream>
#include <random>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

cds::Bounds cds::angleBounds(const DihedralAngleData& metadata)
{
    double defaultAngle = metadata.default_angle_value_;
    return {defaultAngle - metadata.lower_deviation_, defaultAngle + metadata.upper_deviation_};
}

DihedralAngleDataVector cds::likelyMetadata(const DihedralAngleDataVector& entries)
{
    DihedralAngleDataVector returningMetadata;
    returningMetadata.reserve(entries.size());
    for (auto& entry : entries)
    {
        if (entry.weight_ >= 0.01) // HARDCODE EVERYTHING.
        {
            returningMetadata.push_back(entry);
        }
    }
    return returningMetadata;
}

std::string cds::likelyName(const DihedralAngleDataVector& entries)
{
    auto likely = likelyMetadata(entries);
    return likely.empty() ? "Boo" : likely[0].dihedral_angle_name_;
}

cds::AngleWithMetadata cds::defaultAngle(const DihedralAngleData& entry)
{
    return {entry.default_angle_value_, entry};
}

cds::AngleWithMetadata cds::randomDihedralAngleWithinMetadataRange(const DihedralAngleData& entry)
{
    Bounds bounds = angleBounds(entry);
    std::uniform_real_distribution<> angle_distribution(bounds.lower, bounds.upper); // define the range
    double random_angle = angle_distribution(rng);
    return {random_angle, entry};
}

cds::AngleWithMetadata cds::randomAngleEntryUsingMetadata(const DihedralAngleDataVector& entries)
{
    // first randomly pick one of the meta data entries. If there is only one entry, it will randomly always be it.
    std::uniform_int_distribution<> distr(0, (entries.size() - 1)); // define the range
    DihedralAngleData entry = entries.at(distr(rng));
    return randomDihedralAngleWithinMetadataRange(entry);
}
