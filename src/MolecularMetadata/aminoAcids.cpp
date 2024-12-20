#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <string>
#include <vector>
#include <stdexcept>

namespace
{
    using MolecularMetadata::AminoAcid;
    using MolecularMetadata::MoleculeDefinition;

    typedef std::vector<std::pair<std::string, std::string>> BondVector;
    const BondVector backboneBonds {
        { "N", "CA"},
        {"CA",  "C"},
        { "C",  "O"}
    };
    const BondVector carboxylBonds {
        {"C", "OXT"}
    };
    const BondVector carboxylBackboneBonds = codeUtils::vectorAppend(backboneBonds, carboxylBonds);

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

    MoleculeDefinition definitionWithBackbone(const BondVector& backboneBonds, const BondVector& sidechainBonds)
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
            result.push_back({definitionWithBackbone(backboneBonds, sidechain.second),
                              definitionWithBackbone(carboxylBackboneBonds, sidechain.second)});
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
    // clang-format on

    const std::vector<std::string> names     = namesOnly(sidechainBonds);
    const std::vector<AminoAcid> definitions = withBackbones(sidechainBonds);

} // namespace

const std::vector<std::string>& MolecularMetadata::aminoAcidNames()
{
    return names;
}

const std::vector<AminoAcid>& MolecularMetadata::aminoAcids()
{
    return definitions;
}

const AminoAcid& MolecularMetadata::aminoAcid(const std::string& name)
{
    size_t index = codeUtils::indexOf(aminoAcidNames(), name);
    if (index >= aminoAcidNames().size())
    {
        throw std::runtime_error("Oliver! Each query residue residue should be in biology::proteinResidueNames (you "
                                 "need to confirm this before calling this code), or if it's there then it needs an "
                                 "entry in this table. Git gud son. This was the problem residue: " +
                                 name);
    }
    return definitions[index];
}
