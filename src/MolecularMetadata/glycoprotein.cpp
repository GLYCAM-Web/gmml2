#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <unordered_map>

namespace
{
    using namespace glycoproteinMetadata;

    struct GlycosylationEntry
    {
        std::string residue;
        std::vector<std::string> atoms;
        std::vector<SuperimpositionValues> values;
    };

    // Dear future self, the order that you add the atoms to the residue matters for superimposition
    // ie N, CA, CB, not CB, CA, N.
    const std::vector<GlycosylationEntry> glycosylationTableData {
        {"ASN", {"ND2", "CG", "OD1"}, {{109.3, 180, 1.53}, {109.3, 261, 1.325}, {126, 0, 1.22}}},
        {"THR",  {"OG1", "CB", "CA"},  {{112, 68, 1.46}, {109.3, 75, 1.53}, {109.3, 125, 1.53}}},
        {"SER",   {"OG", "CB", "CA"},  {{112, 68, 1.46}, {109.3, 75, 1.53}, {109.3, 125, 1.53}}},
        {"TYR",  {"OH", "CZ", "CE1"},      {{112, 68, 1.46}, {117, 60, 1.35}, {120, 180, 1.37}}}
    };

    std::function<std::string(const GlycosylationEntry&)> getGlycosylationResidue = [](const GlycosylationEntry& a)
    {
        return a.residue;
    };
    std::function<std::vector<std::string>(const GlycosylationEntry&)> getGlycosylationAtoms =
        [](const GlycosylationEntry& a)
    {
        return a.atoms;
    };
    std::function<std::vector<SuperimpositionValues>(const GlycosylationEntry&)> getGlycosylationValues =
        [](const GlycosylationEntry& a)
    {
        return a.values;
    };

    const GlycosylationTable glycosylationTable {codeUtils::vectorMap(getGlycosylationResidue, glycosylationTableData),
                                                 codeUtils::vectorMap(getGlycosylationAtoms, glycosylationTableData),
                                                 codeUtils::vectorMap(getGlycosylationValues, glycosylationTableData)};
} // namespace

std::string glycoproteinMetadata::LookupCodeForAminoAcidName(const std::string queryName)
{
    static const std::unordered_map<std::string, std::string> aminoAcidNameToCodeMap({
        {"ALA", "A"},
        {"ARG", "R"},
        {"ASN", "N"},
        {"ASP", "D"},
        {"CYS", "C"},
        {"GLN", "Q"},
        {"GLU", "E"},
        {"GLY", "G"},
        {"HIS", "H"},
        {"ILE", "I"},
        {"LEU", "L"},
        {"LYS", "K"},
        {"MET", "M"},
        {"PHE", "F"},
        {"PRO", "P"},
        {"SER", "S"},
        {"THR", "T"},
        {"TRP", "W"},
        {"TYR", "Y"},
        {"VAL", "V"},
        {"NLN", "N"},
        {"OLT", "T"},
        {"OLS", "S"},
        {"OLY", "Y"},
        {"CLW", "W"}
    });
    return codeUtils::FindStringInStringMap(queryName, aminoAcidNameToCodeMap);
}

std::string glycoproteinMetadata::LookupLinkTypeForAminoAcidName(const std::string queryName)
{
    static const std::unordered_map<std::string, std::string> residueLinkMap({
        {"ASN", "nLink"},
        {"THR", "oLink"},
        {"SER", "oLink"},
        {"TYR", "oLink"},
        {"TRP", "cLink"},
        {"NLN", "nLink"},
        {"OLT", "oLink"},
        {"OLS", "oLink"},
        {"OLY", "oLink"},
        {"CLW", "cLink"}
    });
    return codeUtils::FindStringInStringMap(queryName, residueLinkMap);
}

std::string glycoproteinMetadata::ConvertGlycosylatedResidueName(const std::string queryname)
{
    static const std::unordered_map<std::string, std::string> GlycamGlycosylatedResidueNameMap({
        {"ASN", "NLN"},
        {"SER", "OLS"},
        {"THR", "OLT"},
        {"TYR", "OLY"},
        {"NLN", "ASN"},
        {"OLS", "SER"},
        {"OLT", "THR"},
        {"OLY", "TYR"}
    });
    return codeUtils::FindStringInStringMap(queryname, GlycamGlycosylatedResidueNameMap);
}

std::string glycoproteinMetadata::GetGlycositeConnectionAtomName(const std::string queryname)
{
    static const std::unordered_map<std::string, std::string> connectionAtomNameMap({
        {"NLN", "ND2"},
        {"ASN", "ND2"},
        {"OLT", "OG1"},
        {"THR", "OG1"},
        {"OLS",  "OG"},
        {"SER",  "OG"},
        {"OLY",  "OH"},
        {"TYR",  "OH"}
    });
    return codeUtils::FindStringInStringMap(queryname, connectionAtomNameMap);
}

const GlycosylationTable& glycoproteinMetadata::defaultGlycosylationTable()
{
    return glycosylationTable;
}
