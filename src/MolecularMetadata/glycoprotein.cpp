#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/find.hpp"

#include <unordered_map>

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
