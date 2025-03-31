#ifndef INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCOPROTEIN_HPP

#include <string>
#include <vector>

namespace glycoproteinMetadata
{
    struct SuperimpositionValues
    {
        double angle;
        double dihedral;
        double distance;
    };

    struct GlycosylationTable
    {
        std::vector<std::string> residueNames;
        std::vector<std::vector<std::string>> atomNames;
        std::vector<std::vector<SuperimpositionValues>> values;
    };

    std::string LookupCodeForAminoAcidName(const std::string queryName);
    std::string LookupLinkTypeForAminoAcidName(const std::string queryName);
    std::string ConvertGlycosylatedResidueName(const std::string queryname);
    std::string GetGlycositeConnectionAtomName(const std::string queryname);
    const GlycosylationTable& defaultGlycosylationTable();
} // namespace glycoproteinMetadata

#endif
