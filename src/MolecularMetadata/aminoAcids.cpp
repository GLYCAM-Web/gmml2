#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <array>
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

    std::vector<std::vector<std::array<std::string, 4>>>
    dihedralAtomsInOrder(const std::vector<std::string>& names, const std::vector<AminoAcid>& definitions,
                         const std::vector<std::pair<std::string, std::vector<std::string>>>& entries)
    {
        std::vector<std::vector<std::array<std::string, 4>>> result;
        result.resize(names.size(), {});
        for (auto& entry : entries)
        {
            const std::string& name = entry.first;
            size_t index            = codeUtils::indexOf(names, name);
            if (index >= names.size())
            {
                throw std::runtime_error("unknown amino acid: " + name);
            }
            const AminoAcid& definition = definitions[index];
            result[index].reserve(entry.second.size());
            for (auto& dihedral : entry.second)
            {
                std::vector<std::string> split = codeUtils::split(dihedral, '-');
                if (split.size() != 4)
                {
                    throw std::runtime_error("error in amino acid metadata for " + name + ", dihedral needs 4 atoms");
                }
                std::array<std::string, 4> atoms {split[0], split[1], split[2], split[3]};
                for (auto& atom : atoms)
                {
                    if (!codeUtils::contains(definition.standard.names, atom))
                    {
                        throw std::runtime_error("unknown atom " + atom + " in dihedral metadata for amino acid " +
                                                 name);
                    }
                }
                result[index].push_back(atoms);
            }
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

    const std::vector<std::string> names     = namesOnly(sidechainBonds);
    const std::vector<AminoAcid> definitions = withBackbones(sidechainBonds);
    const std::vector<std::vector<std::array<std::string, 4>>> dihedralAtoms =
        dihedralAtomsInOrder(names, definitions, sidechainDihedralAtoms);

} // namespace

const std::vector<std::string>& MolecularMetadata::aminoAcidNames()
{
    return names;
}

const std::vector<AminoAcid>& MolecularMetadata::aminoAcids()
{
    return definitions;
}

const AminoAcid& MolecularMetadata::aminoAcid(size_t index)
{
    return definitions[index];
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
    return aminoAcid(index);
}

const std::vector<std::array<std::string, 4>>& MolecularMetadata::aminoAcidDihedrals(size_t index)
{
    return dihedralAtoms[index];
}
