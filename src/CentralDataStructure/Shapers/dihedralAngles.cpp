#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <sstream>

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
