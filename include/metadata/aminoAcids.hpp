#ifndef INCLUDES_METADATA_AMINOACIDS_HPP
#define INCLUDES_METADATA_AMINOACIDS_HPP

#include "include/metadata/elements.hpp"

#include <array>
#include <string>
#include <vector>

namespace gmml
{
    static const std::vector<std::string> proteinResidueNames = {
        "ALA",  "ASP",  "ASN",  "ARG",  "GLY",  "GLU",   "GLN",   "PRO",  "HIS",  "HIP",  "CYS",  "VAL",
        "LEU",  "THR",  "SER",  "LYS",  "MET",  "MSE",   "TYR",   "TRP",  "PHE",  "SEC",  "ILE",  "CYX",
        "CYM",  "HID",  "HIE",  "NLN",  "OLY",  "OLS",   "OLT",   "ARN",  "ASH",  "GLH",  "HYP",  "LYN",
        "NALA", "NASP", "NASN", "NARG", "NGLY", "NGLU",  "NGLN",  "NPRO", "NHIS", "NCYS", "NVAL", "NLEU",
        "NTHR", "NSER", "NLYS", "NMET", "NTYR", "NTRP",  "NPHE",  "NSEC", "NILE", "NCYX", "NCYM", "NHID",
        "NHIE", "NASH", "NGLH", "NHYP", "NLYN", "CNALA", "CNASP", "CASN", "CARG", "CGLY", "CGLU", "CGLN",
        "CPRO", "CHIS", "CCYS", "CVAL", "CLEU", "CTHR",  "CSER",  "CLYS", "CMET", "CTYR", "CTRP", "CPHE",
        "CSEC", "CILE", "CCYX", "CCYM", "CHID", "CHIE",  "CASH",  "CGLH", "CHYP", "CLYN"};

    typedef std::vector<std::pair<std::string, std::string>> BondVector;

    struct AminoAcidTable
    {
        std::vector<std::string> names;
        std::vector<std::string> originalName;
        std::vector<ChemicalFormula> formulas;
        std::vector<std::vector<std::string>> atomNames;
        std::vector<BondVector> bonds;
        std::vector<std::vector<std::array<std::string, 4>>> sidechainDihedralAtoms;
    };

    AminoAcidTable aminoAcidTable();
    ChemicalFormula aminoAcidTerminalAtoms();
    size_t aminoAcidIndex(const AminoAcidTable& table, const std::string& name);
    std::string originalResidueName(const std::string& str);
    const BondVector& carboxylBonds();
} // namespace gmml

#endif
