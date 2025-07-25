#ifndef INCLUDE_METADATA_GLYCAM06RESIDUENAMEGENERATOR_HPP
#define INCLUDE_METADATA_GLYCAM06RESIDUENAMEGENERATOR_HPP

#include <map>
#include <string>
#include <vector>

namespace gmml
{
    namespace metadata
    {
        std::string Glycam06ResidueNameGenerator(
            const std::string& linkages,
            const std::string& preIsomerModifier,
            const std::string& isomer,
            const std::string& inputResName,
            const std::string& ringType,
            const std::string& residueModifier,
            const std::string& configuration);
    } // namespace metadata
} // namespace gmml

#endif
