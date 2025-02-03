#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUETYPES_HPP

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
} // namespace cds
#endif
