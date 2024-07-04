#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06PREPTOSUGARDETAIL_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06PREPTOSUGARDETAIL_HPP

#include <string>
#include <map>
#include <vector>

namespace GlycamMetadata
{
    struct residueMetadata
    {
        std::string isomer          = "";
        std::string resname         = "";
        std::string ringType        = "";
        std::string residueModifier = "";
        std::string configuration   = "";

        std::string toString()
        {
            return isomer + "_" + resname + "_" + ringType + "_" + residueModifier + "_" + configuration;
        }
    };

    residueMetadata Glycam06PrepNameToDetails(const std::string& prepName);
} // namespace GlycamMetadata
#endif