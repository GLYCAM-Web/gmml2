#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"

#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <sstream>
#include <cmath>

using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

std::vector<double> cds::evenlySpaced(double lower, double upper, double approximateIncrement)
{
    double range     = upper - lower;
    int steps        = std::ceil(range / approximateIncrement);
    double increment = range / steps;
    std::vector<double> result;
    result.reserve(steps + 1);
    for (int k = 0; k < steps + 1; k++)
    {
        result.push_back(lower + k * increment);
    }
    return result;
}

std::vector<double> cds::evenlySpacedAngles(double deviation, double increment, const DihedralAngleData& metadata)
{
    double def = metadata.default_angle_value_;
    return cds::evenlySpaced(def - deviation * metadata.lower_deviation_, def + deviation * metadata.upper_deviation_,
                             increment);
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
