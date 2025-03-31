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
        std::vector<std::string> renamedResidues;
        std::vector<std::string> connectingAtomNames;
        std::vector<std::vector<std::string>> atomNames;
        std::vector<std::vector<SuperimpositionValues>> values;
    };

    struct AminoAcidLinkTable
    {
        std::vector<std::string> names;
        std::vector<std::string> codes;
        std::vector<std::string> linkTypes;
    };

    const GlycosylationTable& defaultGlycosylationTable();
    const AminoAcidLinkTable& defaultAminoAcidLinkTable();
} // namespace glycoproteinMetadata

#endif
