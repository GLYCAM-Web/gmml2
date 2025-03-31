#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <unordered_map>

namespace
{
    using namespace glycoproteinMetadata;

    struct GlycosylationEntry
    {
        std::string residue;
        std::string renamed;
        std::string connectingAtom;
        std::vector<std::string> atoms;
        std::vector<SuperimpositionValues> values;
    };

    // Dear future self, the order that you add the atoms to the residue matters for superimposition
    // ie N, CA, CB, not CB, CA, N.
    const std::vector<GlycosylationEntry> glycosylationTableData {
        {"ASN", "NLN", "ND2", {"ND2", "CG", "OD1"}, {{109.3, 180, 1.53}, {109.3, 261, 1.325}, {126, 0, 1.22}}},
        {"THR", "OLT", "OG1",  {"OG1", "CB", "CA"},  {{112, 68, 1.46}, {109.3, 75, 1.53}, {109.3, 125, 1.53}}},
        {"SER", "OLS",  "OG",   {"OG", "CB", "CA"},  {{112, 68, 1.46}, {109.3, 75, 1.53}, {109.3, 125, 1.53}}},
        {"TYR", "OLY",  "OH",  {"OH", "CZ", "CE1"},      {{112, 68, 1.46}, {117, 60, 1.35}, {120, 180, 1.37}}}
    };

    std::function<std::string(const GlycosylationEntry&)> getGlycosylationResidue = [](const GlycosylationEntry& a)
    {
        return a.residue;
    };
    std::function<std::string(const GlycosylationEntry&)> getGlycosylationRenamedResidue =
        [](const GlycosylationEntry& a)
    {
        return a.renamed;
    };
    std::function<std::string(const GlycosylationEntry&)> getGlycosylationConnectingAtom =
        [](const GlycosylationEntry& a)
    {
        return a.connectingAtom;
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

    const GlycosylationTable glycosylationTable {
        codeUtils::vectorMap(getGlycosylationResidue, glycosylationTableData),
        codeUtils::vectorMap(getGlycosylationRenamedResidue, glycosylationTableData),
        codeUtils::vectorMap(getGlycosylationConnectingAtom, glycosylationTableData),
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

const GlycosylationTable& glycoproteinMetadata::defaultGlycosylationTable()
{
    return glycosylationTable;
}
