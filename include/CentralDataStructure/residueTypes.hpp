#ifndef INCLUDE_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP

#include <array>
#include <string>

namespace gmml
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

    struct ResidueAttributes
    {
        ResidueType type = ResidueType::Undefined;
        std::string name = "";          // Gal, Glc,
        std::string glycamCode = "";    // 0GA, OME,
        std::string linkage = "";       // 1-4
        std::string ringType = "";      // p/f
        std::string configuration = ""; // a/b
        std::string isomer = "";        // D/L
        std::string preIsomerModifier = "";
        std::string modifier = "";
        bool isInternal = false; // only care for 2-8. Whether one of the residues has no children.
    };
} // namespace gmml

#endif
