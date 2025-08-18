#ifndef INCLUDE_METADATA_GLYCOPROTEIN_HPP
#define INCLUDE_METADATA_GLYCOPROTEIN_HPP

#include <string>
#include <vector>

namespace gmml
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

    GlycosylationTable defaultGlycosylationTable();
    AminoAcidLinkTable defaultAminoAcidLinkTable();
} // namespace gmml

#endif
