#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP

#include <string>
#include <map>
#include <vector>

namespace GlycamMetadata
{
    std::string Glycam06ResidueNameGenerator(std::string linkages, std::string isomer, std::string inputResName,
                                             std::string ringType, std::string residueModifier,
                                             std::string configuration);
} // namespace GlycamMetadata
#endif
