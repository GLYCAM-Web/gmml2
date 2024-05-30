#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP

#include <string>
#include <vector>

namespace glycoproteinMetadata
{
    std::string LookupCodeForAminoAcidName(const std::string queryName);
    std::string LookupLinkTypeForAminoAcidName(const std::string queryName);
    std::string ConvertGlycosylatedResidueName(const std::string queryname);
    std::string GetGlycositeConnectionAtomName(const std::string queryname);
} // namespace glycoproteinMetadata

#endif
