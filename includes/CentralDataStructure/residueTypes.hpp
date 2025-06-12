#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP

#include <string>
#include <array>

namespace cds
{
    enum ResidueType
    {
        Protein,
        Sugar,
        Aglycone,
        Derivative,
        Solvent,
        Deoxy,
        ProteinCappingGroup,
        Undefined,
        ResidueTypeCount
    };

    // For dihedralangledata.cpp
    inline std::string residueTypeToString(cds::ResidueType resType)
    {
        const std::array<std::string, 9> residueTypeStrings = {
            "amino-acid", "monosaccharide",      "aglycon",   "derivative",      "solvent",
            "deoxy",      "proteinCappingGroup", "undefined", "ResidueTypeCount"};
        size_t index = static_cast<size_t>(resType);
        return (index < residueTypeStrings.size()) ? residueTypeStrings[index] : "UNKNOWN";
    };

    struct ResidueAttributes
    {
        ResidueType type              = ResidueType::Undefined;
        std::string name              = ""; // Gal, Glc,
        std::string linkage           = ""; // 1-4
        std::string ringType          = ""; // p/f
        std::string configuration     = ""; // a/b
        std::string isomer            = ""; // D/L
        std::string preIsomerModifier = "";
        std::string modifier          = "";
        std::string glycamCode        = "";    // 0GA, OME,
        bool isInternal               = false; // only care for 2-8. Whether one of the residues has no children.
    };
} // namespace cds
#endif
