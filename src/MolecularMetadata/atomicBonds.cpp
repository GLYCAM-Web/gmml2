#include "includes/MolecularMetadata/atomicBonds.hpp"

#include <array>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace
{
    using MolecularMetadata::Element;
    typedef std::map<std::array<Element, 2>, std::pair<double, double>> BondLengthMap;

    const double maxCutOff = 1.65;
    const double minCutOff = 0.7;
    const double defaultBondLength = 1.4;
    const double hydrogenMaxCutoff = 1.15;

    BondLengthMap bondLengthRanges()
    {
        Element C = Element::C;
        Element O = Element::O;
        Element N = Element::N;
        Element P = Element::P;
        Element S = Element::S;
        return {
            {{C, C},  {1.22, 1.67}},
            {{C, O},  {1.08, 1.68}},
            {{O, C},  {1.08, 1.68}},
            {{C, N},  {1.26, 1.55}},
            {{N, C},  {1.26, 1.55}},
            {{O, P}, {1.35, 1.776}},
            {{P, O}, {1.35, 1.776}},
            {{O, S},  {1.43, 1.78}},
            {{S, O},  {1.43, 1.78}},
            {{N, S},  {1.62, 1.77}},
            {{S, N},  {1.62, 1.77}},
            {{C, S},  {1.50, 1.91}},
            {{S, C},  {1.50, 1.91}},
            {{S, S},  {1.50, 2.10}}
        };
    }

    const BondLengthMap bondLengthRangeMap = bondLengthRanges();

    std::pair<double, double> bondLengthRange(Element atom1Element, Element atom2Element)
    { // Using PDB bond length statistics provided by Chenghua on 2/5/19

        std::array<Element, 2> bothAtoms = {atom1Element, atom2Element};

        if (bondLengthRangeMap.find(bothAtoms) != bondLengthRangeMap.end())
        {
            return bondLengthRangeMap.at(bothAtoms);
        }
        else
        {
            // gmml::log(__LINE__, __FILE__,  gmml::INF, "Using default binding cutoff of 1.65");
            return {minCutOff, maxCutOff};
        }
    }

    struct SpecificBondLength
    {
        std::string type1_; // One of the atom types
        std::string type2_; // The other atom type
        double length_;     // Default length in Angstroms
    };

    // clang-format off
    const std::vector<SpecificBondLength> specificBondLengths = {
        {"Cg",  "N",  1.450}, // Copy of Cg-Ng from GLYCAM06
        {"CT", "Os",  1.410}, // Copy of CT-OS from parm10.dat
        {"2C", "Os",  1.410}, // Copy of CT-OS from parm10.dat
        {"3C", "Os",  1.410}, // Copy of CT-OS from parm10.dat
        { "S", "Ng",  1.675}, // N-Sulfate - Using avg value from ZULPIF and ZULPIF01 (CSD 1.638 A)
        {"Cg", "Sm",  1.810}, // Changed from 222.0 based on methanethiol nmodes
        {"N3",  "H",  1.010}, // Parm10
        {"N3", "Cg",  1.490}, // Methanaminium (eqm value from crystal avg)
        {"Cj", "Ck",  1.467}, // CRC manual for 1,3-butadiene
        {"Os", "Cj",  1.359}, // copy of Os-Ck
        {"Os", "Ck",  1.359}, // Methoxyethene (JACS 1993, 115, 11921)
        {"Cj", "Cj",  1.337}, // copy of Ck-Ck
        {"Ck", "Ck",  1.337}, // JCC 1996, 17 (5&6),669
        {"Cg", "Cj",  1.514}, // copy of Cg-Ck
        {"Cg", "Ck",  1.514}, // JCC 1996, 17 (5&6),669
        {"Cj", "Ha",  1.095}, // copy of Ck-Ha
        {"Ck", "Ha",  1.095}, // Ethane for alkenes
        {"Cg", "Hp",  1.095}, // Copy of Cg-Hc
        {"NT",  "H",  1.010}, // Parm99
        {"NT", "Cg",  1.470}, // K calculated from methyl amine (eqm value from crystal average)
        {"Cg", "Cp",  1.520}, // Copy of Cg-Cg
        {"Cp", "H2",  1.090}, // Copy of Cg-H1
        {"Cp", "H1",  1.090}, // Copy of Cg-H1
        {"Cp", "Os",  1.460}, // Copy of Cg-Os
        { "P", "Os",  1.610}, // Parm94
        { "P", "O2",  1.480}, // Parm10
        { "S", "Os",  1.589}, // K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)
        { "S", "O2",  1.440}, // K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)
        { "C", "Os",  1.323}, // Parm99
        {"OW", "HW", 0.9572}, // TIP3P water
        {"HW", "HW", 1.5136}, // TIP3P water
        {"Cg", "Cg",  1.520}, // Butane (gauche, and trans)
        {"Cg", "Hc",  1.090}, // Parm94
        {"Cg", "H1",  1.090}, // Parm94
        {"Cg", "H2",  1.090}, // Parm94
        {"Cg", "Oh",  1.430}, // Methanol
        {"Cg",  "C",  1.530}, // 2-Methylpropanoate
        {"Cg", "Ng",  1.450}, // N-Methylethanamide
        { "C",  "O",  1.229}, // Parm10
        { "C", "Ng",  1.335}, // Parm94
        { "C", "O2",  1.250}, // Parm10
        { "C", "Hc",  1.090}, // Parm91
        { "C", "H1",  1.092}, // Methanol
        {"Oh", "Ho",  0.960}, // Methanol
        {"Ng",  "H",  1.010}, // Parm94
        {"Cy", "Oh",  1.410}, // Parm94 - for sialic acid only!
        {"Cg", "Os",  1.460}, // Parm94   K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)
        {"Cg", "Oy",  1.410}, // Parm94 - for sialic acid only!
        {"Cy", "Cg",  1.520}, // Butane (gauche, and trans) - for sialic acid only!
        {"Cy", "Os",  1.410}, // Parm94 - for sialic acid only!
        {"Cy", "Oy",  1.410}, // Parm94 - for sialic acid only!
        {"Cy",  "C",  1.530}  // 2-Methylpropanoate - for sialic acid only!
    };
    // clang-format on
} // namespace

double MolecularMetadata::hydrogenCovalentBondMaxLength() { return hydrogenMaxCutoff; }

double MolecularMetadata::maxBondLengthByAtomType(Element atom1Element, Element atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::pair<double, double> result = bondLengthRange(atom1Element, atom2Element);
    return result.second;
}

double MolecularMetadata::specificBondLength(const std::string& query1, const std::string& query2)
{
    for (auto& entry : specificBondLengths)
    { // Search bidirectionally e.g Cg-Os, Os-Cg
        if (((entry.type1_ == query1) && (entry.type2_ == query2)) ||
            ((entry.type1_ == query2) && (entry.type2_ == query1)))
        {
            return entry.length_;
        }
    }
    return defaultBondLength;
}
