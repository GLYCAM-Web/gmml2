#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06FUNCTIONS_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06FUNCTIONS_HPP
#include <map>
#include <string>
#include <vector>

namespace GlycamMetadata
{
    std::string GetGlycam06ResidueLinkageCode(const std::string& query);
    std::string GetNameForCode(const std::string& query);
    std::string GetCodeForName(const std::string& query);
    std::string GetTypeForCode(const std::string& query);
    std::string GetDescriptiveNameForGlycamResidueName(const std::string& residueNameInGLYCAMFormat);
    double GetAdjustmentCharge(const std::string& queryResidueCode);
    std::string GetAdjustmentAtom(const std::string& queryResidueCode);
    std::string GetConnectionAtomForResidue(const std::string& query);
} // namespace GlycamMetadata
#endif
