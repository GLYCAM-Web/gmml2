#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP

#include <string>
#include <map>
#include <vector>

namespace GlycamMetadata
{
    std::string Glycam06ResidueNameGenerator(const std::string& linkages, const std::string& isomer,
                                             const std::string& inputResName, const std::string& ringType,
                                             const std::string& residueModifier, const std::string& configuration);
} // namespace GlycamMetadata
#endif
