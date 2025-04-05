#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <array>
#include <vector>
#include <stdexcept>

namespace
{
    using MolecularMetadata::BondVector;

    struct AminoAcid
    {
        std::vector<std::string> atomNames;
        BondVector bonds;
    };

    const std::vector<std::pair<std::string, double>> molecularWeight {
        {"ALA",  89},
        {"ARG", 174},
        {"ASN", 132},
        {"ASP", 133},
        {"CYS", 121},
        {"CYX", 120},
        {"GLN", 146},
        {"GLU", 147},
        {"GLY",  75},
        {"HIS", 155},
        {"ILE", 131},
        {"LEU", 131},
        {"LYS", 146},
        {"MET", 149},
        {"MSE",   0},
        {"PHE", 165},
        {"PRO", 115},
        {"SEC",   0},
        {"SER", 105},
        {"THR", 119},
        {"TRP", 204},
        {"TYR", 181},
        {"VAL", 117}
    };

    const std::vector<std::pair<std::string, std::string>> originalResidueNames = {
        {"HIE", "HIS"},
        {"HID", "HIS"},
        {"HIP", "HIS"},
        {"CYX", "CYS"},
        {"CYM", "CYS"},
        {"NLN", "ASN"},
        {"OLS", "SER"},
        {"OLT", "THR"},
        {"OLY", "TYR"},
        {"ARN", "ARG"},
        {"ASH", "ASP"},
        {"GLH", "GLU"},
        {"LYN", "LYS"},
    };

    const BondVector backboneBondVector {
        { "N", "CA"},
        {"CA",  "C"},
        { "C",  "O"}
    };
    const BondVector carboxylBondVector {
        {"C", "OXT"}
    };

    std::vector<std::string> uniqueAtomNamesPresentInBonds(const BondVector& bonds)
    {
        std::vector<std::string> atomNames;
        atomNames.reserve(bonds.size() * 2);
        for (auto& bond : bonds)
        {
            atomNames.push_back(bond.first);
            atomNames.push_back(bond.second);
        }
        return codeUtils::sorted(codeUtils::uniqueOnly(atomNames));
    }

    AminoAcid definitionWithBackbone(const BondVector& backboneBonds, const BondVector& sidechainBonds)
    {
        BondVector bonds = codeUtils::vectorAppend(backboneBonds, sidechainBonds);
        return {uniqueAtomNamesPresentInBonds(bonds), bonds};
    }

    std::vector<AminoAcid> withBackbones(const std::vector<std::pair<std::string, BondVector>>& sidechains)
    {
        std::vector<AminoAcid> result;
        result.reserve(sidechains.size());
        for (auto& sidechain : sidechains)
        {
            result.push_back(definitionWithBackbone(backboneBondVector, sidechain.second));
        }
        return result;
    }

    std::vector<std::string> namesOnly(const std::vector<std::pair<std::string, BondVector>>& sidechains)
    {
        std::vector<std::string> result;
        result.reserve(sidechains.size());
        for (auto& sidechain : sidechains)
        {
            result.push_back(sidechain.first);
        }
        return result;
    }

    std::vector<std::array<std::string, 4>> splitDihedralAtoms(const std::string& name,
                                                               const std::vector<std::string>& dihedrals)
    {
        std::vector<std::array<std::string, 4>> result;
        result.reserve(dihedrals.size());
        for (auto& dihedral : dihedrals)
        {
            std::vector<std::string> split = codeUtils::split(dihedral, '-');
            if (split.size() != 4)
            {
                throw std::runtime_error("error in amino acid metadata for " + name + ", dihedral needs 4 atoms");
            }
            result.push_back({split[0], split[1], split[2], split[3]});
        }
        return result;
    }

    // clang-format off
    const std::vector<std::pair<std::string, BondVector>> sidechainBonds {
        {"ALA", {{"CA", "CB"}}},
        {"ARN", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "NE"}, {"NE", "CZ"}, {"CZ", "NH1"}, {"CZ", "NH2"}}},
        {"ARG", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "NE"}, {"NE", "CZ"}, {"CZ", "NH1"}, {"CZ", "NH2"}}},
        {"ASH", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "OD2"}}},
        {"ASP", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "OD2"}}},
        {"ASN", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "ND2"}}},
        {"NLN", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "ND2"}}},
        {"CYM", {{"CA", "CB"}, {"CB", "SG"}}},
        {"CYS", {{"CA", "CB"}, {"CB", "SG"}}},
        {"CYX", {{"CA", "CB"}, {"CB", "SG"}}},
        {"GLH", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "OE2"}}},
        {"GLU", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "OE2"}}},
        {"GLN", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "NE2"}}},
        {"GLY", {}},
        {"HID", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD2"}, {"CD2", "NE2"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}}},
        {"HIE", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD2"}, {"CD2", "NE2"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}}},
        {"HIP", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD2"}, {"CD2", "NE2"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}}},
        {"HIS", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD2"}, {"CD2", "NE2"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}}},
        {"ILE", {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}, {"CG1", "CD1"}}},
        {"LEU", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}}},
        {"LYN", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}, {"CE", "NZ"}}},
        {"LYS", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}, {"CE", "NZ"}}},
        {"MET", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "SD"}, {"SD", "CE"}}},
        {"PHE", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "CE1"}, {"CE1", "CZ"}, {"CG", "CD2"}, {"CD2", "CE2"}, {"CE2", "CZ"}}},
        {"PRO", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "N"}}},
        {"OLS", {{"CA", "CB"}, {"CB", "OG"}}},
        {"SER", {{"CA", "CB"}, {"CB", "OG"}}},
        {"OLT", {{"CA", "CB"}, {"CB", "OG1"}, {"CB", "CG2"}}},
        {"THR", {{"CA", "CB"}, {"CB", "OG1"}, {"CB", "CG2"}}},
        {"TRP", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}, {"CD1", "NE1"}, {"CD2", "CE2"}, {"CD2", "CE3"}, {"NE1", "CE2"}, {"CE2", "CZ2"}, {"CE3", "CZ3"}, {"CZ2", "CH2"}, {"CZ3", "CH2"}}},
        {"OLY", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "CE1"}, {"CE1", "CZ"}, {"CG", "CD2"}, {"CD2", "CE2"}, {"CE2", "CZ"}, {"CZ", "OH"}}},
        {"TYR", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "CE1"}, {"CE1", "CZ"}, {"CG", "CD2"}, {"CD2", "CE2"}, {"CE2", "CZ"}, {"CZ", "OH"}}},
        {"VAL", {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}}},
        {"SEC", {{"CA", "CB"}, {"CB", "SE"}}},
        {"MSE", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "SE"}, {"SE", "CE"}}}};

    const std::vector<std::pair<std::string, std::vector<std::string>>> sidechainDihedralAtoms {
        {"ARN", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-NE", "CG-CD-NE-CZ", "CD-NE-CZ-NH1"}},
        {"ARG", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-NE", "CG-CD-NE-CZ", "CD-NE-CZ-NH1"}},
        {"ASN", {"N-CA-CB-CG", "CA-CB-CG-OD1"}},
        {"ASH", {"N-CA-CB-CG", "CA-CB-CG-OD1"}},
        {"ASP", {"N-CA-CB-CG", "CA-CB-CG-OD1"}},
        {"CYM", {"N-CA-CB-SG"}},
        {"CYS", {"N-CA-CB-SG"}},
        {"GLH", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-OE1"}},
        {"GLN", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-OE1"}},
        {"GLU", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-OE1"}},
        {"HID", {"N-CA-CB-CG", "CA-CB-CG-ND1"}},
        {"HIE", {"N-CA-CB-CG", "CA-CB-CG-ND1"}},
        {"HIP", {"N-CA-CB-CG", "CA-CB-CG-ND1"}},
        {"HIS", {"N-CA-CB-CG", "CA-CB-CG-ND1"}},
        {"ILE", {"N-CA-CB-CG1", "CA-CB-CG1-CD1"}},
        {"LEU", {"N-CA-CB-CG", "CA-CB-CG-CD1"}},
        {"LYN", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-CE", "CG-CD-CE-NZ"}},
        {"LYS", {"N-CA-CB-CG", "CA-CB-CG-CD", "CB-CG-CD-CE", "CG-CD-CE-NZ"}},
        {"MET", {"N-CA-CB-CG", "CA-CB-CG-SD", "CB-CG-SD-CE"}},
        {"PHE", {"N-CA-CB-CG", "CA-CB-CG-CD1"}},
        {"SER", {"N-CA-CB-OG"}},
        {"THR", {"N-CA-CB-OG1"}},
        {"TRP", {"N-CA-CB-CG", "CA-CB-CG-CD1"}},
        {"TYR", {"N-CA-CB-CG", "CA-CB-CG-CD1"}},
        {"VAL", {"N-CA-CB-CG1"}}
    };

    // clang-format on

    MolecularMetadata::AminoAcidTable createTable()
    {
        MolecularMetadata::AminoAcidTable table;
        size_t count = molecularWeight.size() + originalResidueNames.size();
        table.names.reserve(count);
        table.originalName.reserve(count);
        table.weights.reserve(count);
        for (auto& a : molecularWeight)
        {
            table.names.push_back(a.first);
            table.originalName.push_back(a.first);
            table.weights.push_back(a.second);
        }
        for (auto& a : originalResidueNames)
        {
            table.names.push_back(a.first);
            table.originalName.push_back(a.second);
            size_t originalIndex = codeUtils::indexOf(table.names, a.second);
            table.weights.push_back(table.weights[originalIndex]);
        }
        table.bonds.resize(count);
        table.atomNames.resize(count);
        const std::vector<std::string> atomDefinitionNames = namesOnly(sidechainBonds);
        const std::vector<AminoAcid> atomDefinitions       = withBackbones(sidechainBonds);
        for (size_t n = 0; n < atomDefinitionNames.size(); n++)
        {
            size_t index = codeUtils::indexOf(table.names, atomDefinitionNames[n]);
            if (index == table.names.size())
            {
                throw std::runtime_error(atomDefinitionNames[n] + " from atom definitions not found\n");
            }
            table.atomNames[index] = atomDefinitions[n].atomNames;
            table.bonds[index]     = atomDefinitions[n].bonds;
        }
        table.sidechainDihedralAtoms.resize(count);
        for (auto& a : sidechainDihedralAtoms)
        {
            size_t index = codeUtils::indexOf(table.names, a.first);
            if (index == table.names.size())
            {
                throw std::runtime_error(a.first + " from dihedrals not found\n");
            }
            std::vector<std::array<std::string, 4>> atoms = splitDihedralAtoms(a.first, a.second);
            table.sidechainDihedralAtoms[index]           = atoms;
        }
        return table;
    }

    const MolecularMetadata::AminoAcidTable aminoAcidTableVar = createTable();

} // namespace

const MolecularMetadata::AminoAcidTable& MolecularMetadata::aminoAcidTable()
{
    return aminoAcidTableVar;
}

size_t MolecularMetadata::aminoAcidIndex(const AminoAcidTable& table, const std::string& name)
{
    size_t index = codeUtils::indexOf(table.names, name);
    if (index >= table.names.size())
    {
        throw std::runtime_error("Error: amino acid not found in metadata: " + name);
    }
    return index;
}

std::string MolecularMetadata::originalResidueName(const std::string& str)
{
    auto iter = std::find_if(originalResidueNames.begin(), originalResidueNames.end(),
                             [&](const std::pair<std::string, std::string>& pair)
                             {
                                 return pair.first == str;
                             });
    if (iter != originalResidueNames.end())
    {
        return iter->second;
    }
    return str;
}

const MolecularMetadata::BondVector& MolecularMetadata::carboxylBonds()
{
    return carboxylBondVector;
}
