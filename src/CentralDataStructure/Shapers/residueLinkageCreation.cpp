#include <cmath>
#include <sstream>
#include <iostream>
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/psiAngleHydrogen.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

namespace
{
    cds::ResidueAttributes generateAttributes(const cds::Residue* residue)
    {
        if (residue->GetType() == cds::ResidueType::Protein)
        {
            return cds::ResidueAttributes {cds::ResidueType::Protein, residue->getName(), residue->getName()};
        }
        return residue->getAttributes(); // Only carbs have their attributes set.
    }

    cds::ResidueLinkAttributes toAttributes(const cds::ResidueLink link)
    {
        return {
            {generateAttributes(link.residues.first), generateAttributes(link.residues.second)},
            {            link.atoms.first->getName(),             link.atoms.second->getName()}
        };
    }

    std::string determineLinkageNameFromResidueNames(const cds::ResidueLinkAttributes link)
    {
        std::string residue1Name =
            GlycamMetadata::GetDescriptiveNameForGlycamResidueName(link.residues.first.glycamCode);
        std::string residue2Name =
            GlycamMetadata::GetDescriptiveNameForGlycamResidueName(link.residues.second.glycamCode);
        std::string atom1Name = link.atoms.first;
        std::string atom2Name = link.atoms.second;
        char link1            = *atom1Name.rbegin(); //
        char link2            = *atom2Name.rbegin(); // Messy for Acetyl.
        std::stringstream linkageName;
        linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
        return linkageName.str();
    }

    std::vector<std::vector<size_t>> findResidueLinkageMetadata(cds::ResidueLinkAttributes link)
    {
        std::string firstAtom                 = link.atoms.first;
        std::string secondAtom                = link.atoms.second;
        cds::ResidueAttributes& firstResidue  = link.residues.first;
        cds::ResidueAttributes& secondResidue = link.residues.second;
        std::vector<std::vector<size_t>> matching_entries =
            GlycamMetadata::getDihedralAngleDataEntriesForLinkage(firstAtom, firstResidue, secondAtom, secondResidue);
        if (matching_entries.empty())
        { // Trying the reverse order
            matching_entries = GlycamMetadata::getDihedralAngleDataEntriesForLinkage(secondAtom, secondResidue,
                                                                                     firstAtom, firstResidue);
        }
        if (matching_entries.empty())
        {
            std::stringstream ss;
            ss << "No Metadata entries found for connection between " << firstResidue.glycamCode << "@" << firstAtom
               << " and " << secondResidue.glycamCode << "@" << secondAtom << "\n";
            ss << "Note that order should be reducing atom - anomeric atom, but I've tried reversing the order and it "
                  "didn't fix the issue.\n";
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            throw std::runtime_error(ss.str());
        }
        return matching_entries;
    }

    std::vector<cds::RotatableDihedral>
    createRotatableDihedrals(const std::string& linkageName, const std::vector<cds::DihedralAtoms>& dihedralAtoms,
                             const std::vector<std::vector<size_t>>& metadataIndices)
    {
        std::vector<cds::RotatableDihedral> rotatableDihedrals;
        rotatableDihedrals.reserve(dihedralAtoms.size());
        for (size_t n = 0; n < dihedralAtoms.size(); n++)
        {
            const std::vector<size_t>& currentMetadata = metadataIndices[n];
            if (!currentMetadata.empty())
            {
                rotatableDihedrals.emplace_back(cds::RotatableDihedral {dihedralAtoms[n], {}, 0});
            }
            else
            {
                std::stringstream ss;
                ss << "Problem with the metadata found in gmml for this linkage. No metadata found for dihedral with "
                      "bond: ";
                ss << linkageName << "\n";
                ss << "At dihedral index: " << n << "\n";
                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                throw std::runtime_error(ss.str());
            }
        }
        return rotatableDihedrals;
    }

    void validateRotamerTypes(const GlycamMetadata::DihedralAngleDataTable& table, const cds::ResidueLinkage& linkage)
    {
        for (auto& metadataVector : linkage.dihedralMetadata)
        {
            for (auto& metadata : metadataVector)
            {
                if (table.entries[metadata].rotamer_type_ != linkage.rotamerType)
                {
                    throw std::runtime_error("mismatching rotamer types in residue linkage: " + print(table, linkage));
                }
            }
        }
    }

    void validateConformerMetadata(const GlycamMetadata::DihedralAngleDataTable& table,
                                   const cds::ResidueLinkage& linkage)
    {
        auto& dihedrals       = linkage.rotatableDihedrals;
        auto& metadata        = linkage.dihedralMetadata;
        size_t conformerCount = metadata[0].size();
        for (size_t n = 1; n < dihedrals.size(); n++)
        {
            if (metadata[n].size() != conformerCount)
            {
                throw std::runtime_error("error: different number of conformers in linkage: " +
                                         cds::print(table, linkage));
            }
        }

        for (size_t conformer = 0; conformer < conformerCount; conformer++)
        {
            double weight = table.entries[metadata[0][conformer]].weight_;
            for (size_t n = 1; n < dihedrals.size(); n++)
            {
                if (std::fabs(weight - table.entries[metadata[n][conformer]].weight_) >= 1e-10)
                {
                    throw std::runtime_error("error: conformers have different weight in linkage: " +
                                             cds::print(table, linkage));
                }
            }
        }
    }

    void throwMissingMetadataError(const cds::ResidueLinkAttributes& names,
                                   const std::vector<cds::DihedralAtoms>& dihedralAtoms,
                                   const std::vector<std::vector<size_t>>& metadataIndices)
    {
        size_t dihedralCount = dihedralAtoms.size();
        size_t missingCount  = dihedralCount - metadataIndices.size();
        std::ostringstream stream;
        stream << "Insufficient metadata found for the linkage between " << names.residues.first.name << " and "
               << names.residues.second.name << ". Missing metadata for ";
        if (missingCount > 1)
        {
            stream << "dihedrals " << dihedralCount - missingCount + 1 << "-" << dihedralCount;
        }
        else
        {
            stream << "dihedral " << dihedralCount;
        }
        std::vector<std::string> dihedralInfo;
        for (size_t n = metadataIndices.size(); n < dihedralAtoms.size(); n++)
        {
            std::vector<std::string> atomNames;
            for (auto& atom : dihedralAtoms[n])
            {
                atomNames.push_back(atom->getName());
            }
            dihedralInfo.push_back(codeUtils::join(", ", atomNames));
        }
        stream << ": [" << codeUtils::join("], [", dihedralInfo) << "]";
        throw std::runtime_error(stream.str());
    }
} // namespace

unsigned long long cds::generateResidueLinkageIndex()
{ // static keyword means it is created only once and persists beyond scope of code block.
    static unsigned long long s_ResidueLinkageIndex = 0;
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the
                                    // value in the copy
}

cds::ResidueLink cds::findResidueLink(std::pair<cds::Residue*, cds::Residue*> residues)
{
    std::vector<cds::Atom*> atoms;
    bool found = false;
    cdsSelections::FindAtomsConnectingResidues(residues.first->getAtoms().at(0), residues.first, residues.second,
                                               &atoms, &found);
    if (atoms.size() >= 2)
    {
        return {
            residues, {atoms[0], atoms[1]}
        };
    }
    else
    {
        throw std::runtime_error("Two residues passed into findResidueLink that have no connection atoms.");
    }
}

void cds::determineAtomsThatMove(std::vector<RotatableDihedral>& dihedrals)
{
    for (auto& dihedral : dihedrals)
    {
        auto& atoms = dihedral.atoms;
        std::vector<cds::Atom*> atoms_that_move;
        atoms_that_move.push_back(atoms[2]);
        cdsSelections::FindConnectedAtoms(atoms_that_move, atoms[1]);
        atoms_that_move.erase(atoms_that_move.begin());
        dihedral.movingAtoms = atoms_that_move;
    }
}

cds::ResidueLinkage cds::createResidueLinkage(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                              ResidueLink& link)
{
    std::string firstId                  = cds::residueStringId(link.residues.first);
    std::string secondId                 = cds::residueStringId(link.residues.second);
    ResidueLinkAttributes linkAttributes = toAttributes(link);
    int local_debug                      = 1;
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Maybe Finding connection between " + firstId + " :: " + secondId);
    }
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Finding connection between " + firstId + " :: " + secondId);
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Connection atoms are from: " + link.atoms.first->getId() + " to " + link.atoms.second->getId());
    }
    std::vector<DihedralAtoms> dihedralAtoms = cdsSelections::findRotatableDihedralsConnectingResidues(link);
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Finding metadata for " + firstId + " :: " + secondId);
    }
    std::vector<std::vector<size_t>> metadata = findResidueLinkageMetadata(linkAttributes);
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Metadata found:");
        for (auto& entry : metadata)
        {
            for (auto& dihedralAngleData : entry)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, metadataTable.entries[dihedralAngleData].print());
            }
        }
    }
    if (metadata.size() < dihedralAtoms.size())
    {
        throwMissingMetadataError(linkAttributes, dihedralAtoms, metadata);
    }
    if (metadata.size() > dihedralAtoms.size())
    {
        std::string message =
            "Found metadata for rotatable bonds that do not exist.\nCheck both dihedralangledata metadata "
            "and ResidueLinkage::FindRotatableDihedralsConnectingResidues.\nNote this is normal for a sialic acid "
            "with multiple 2-7, 2-8 and or 2-9 linkages and this warning can be ignored\n";
        gmml::log(__LINE__, __FILE__, gmml::WAR, message);
    }
    // ensure that metadata and dihedral atoms have the same dimensions
    metadata.erase(metadata.begin() + dihedralAtoms.size(), metadata.end());

    auto& residues = link.residues;
    createHydrogenForPsiAngles(metadataTable, residues.second, dihedralAtoms, metadata);
    std::string linkageName                  = determineLinkageNameFromResidueNames(linkAttributes);
    std::vector<RotatableDihedral> dihedrals = createRotatableDihedrals(linkageName, dihedralAtoms, metadata);
    determineAtomsThatMove(dihedrals);

    unsigned long long index = generateResidueLinkageIndex();

    if (dihedrals.empty() || metadata[0].empty())
    {
        throw std::runtime_error("missing dihedrals or metadata in residue linkage: " + print(link));
    }

    const GlycamMetadata::DihedralAngleData& firstMetadata = metadataTable.entries[metadata[0][0]];
    GlycamMetadata::RotamerType rotamerType                = firstMetadata.rotamer_type_;
    const std::vector<std::string>& cond1                  = firstMetadata.residue1_conditions_;
    bool isDerivative                                      = cond1.size() > 0 && cond1[0] == "derivative";
    ResidueLinkage linkage {link, dihedrals, metadata, rotamerType, isDerivative, index, linkageName};

    validateRotamerTypes(metadataTable, linkage);
    if (rotamerType == GlycamMetadata::RotamerType::conformer)
    {
        validateConformerMetadata(metadataTable, linkage);
    }

    return linkage;
}
