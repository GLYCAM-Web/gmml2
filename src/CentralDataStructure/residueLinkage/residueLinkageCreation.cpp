#include "include/CentralDataStructure/residueLinkage/residueLinkageCreation.hpp"

#include "include/CentralDataStructure/Selections/atomSelections.hpp"
#include "include/CentralDataStructure/Selections/residueSelections.hpp"
#include "include/CentralDataStructure/Selections/shaperSelections.hpp"
#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/CentralDataStructure/residueLinkage/psiAngleHydrogen.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageFunctions.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/glycam06Functions.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <cmath>
#include <sstream>

namespace gmml
{
    namespace
    {
        ResidueAttributes generateAttributes(const Residue* residue)
        {
            if (residue->GetType() == ResidueType::Protein)
            {
                return ResidueAttributes {ResidueType::Protein, residue->getName(), residue->getName()};
            }
            return residue->getAttributes(); // Only carbs have their attributes set.
        }

        ResidueLinkAttributes toAttributes(const ResidueLink link)
        {
            return {
                {generateAttributes(link.residues.first), generateAttributes(link.residues.second)},
                {            link.atoms.first->getName(),             link.atoms.second->getName()}
            };
        }

        std::string determineLinkageNameFromResidueNames(const ResidueLinkAttributes link)
        {
            std::string residue1Name = metadata::GetDescriptiveNameForGlycamResidueName(link.residues.first.glycamCode);
            std::string residue2Name =
                metadata::GetDescriptiveNameForGlycamResidueName(link.residues.second.glycamCode);
            std::string atom1Name = link.atoms.first;
            std::string atom2Name = link.atoms.second;
            char link1 = *atom1Name.rbegin(); //
            char link2 = *atom2Name.rbegin(); // Messy for Acetyl.
            std::stringstream linkageName;
            linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
            return linkageName.str();
        }

        std::vector<std::vector<size_t>> findResidueLinkageMetadata(ResidueLinkAttributes link)
        {
            std::string firstAtom = link.atoms.first;
            std::string secondAtom = link.atoms.second;
            ResidueAttributes& firstResidue = link.residues.first;
            ResidueAttributes& secondResidue = link.residues.second;
            std::vector<std::vector<size_t>> matching_entries =
                getDihedralAngleDataEntriesForLinkage(firstAtom, firstResidue, secondAtom, secondResidue);
            if (matching_entries.empty())
            { // Trying the reverse order
                matching_entries =
                    getDihedralAngleDataEntriesForLinkage(secondAtom, secondResidue, firstAtom, firstResidue);
            }
            if (matching_entries.empty())
            {
                std::stringstream ss;
                ss << "No Metadata entries found for connection between " << firstResidue.glycamCode << "@" << firstAtom
                   << " and " << secondResidue.glycamCode << "@" << secondAtom << "\n";
                ss << "Note that order should be reducing atom - anomeric atom, but I've tried reversing the order and "
                      "it "
                      "didn't fix the issue.\n";
                util::log(__LINE__, __FILE__, util::ERR, ss.str());
                throw std::runtime_error(ss.str());
            }
            return matching_entries;
        }

        std::vector<RotatableBond> createRotatableBonds(
            const std::string& linkageName,
            const std::vector<DihedralAtoms>& dihedralAtoms,
            const std::vector<std::vector<size_t>>& metadataIndices)
        {
            std::vector<RotatableBond> result;
            result.reserve(dihedralAtoms.size());
            for (size_t n = 0; n < dihedralAtoms.size(); n++)
            {
                const std::vector<size_t>& currentMetadata = metadataIndices[n];
                if (!currentMetadata.empty())
                {
                    result.emplace_back(RotatableBond {dihedralAtoms[n], {}, 0});
                }
                else
                {
                    std::stringstream ss;
                    ss << "Problem with the metadata found in gmml for this linkage. No metadata found for dihedral "
                          "with "
                          "bond: ";
                    ss << linkageName << "\n";
                    ss << "At dihedral index: " << n << "\n";
                    util::log(__LINE__, __FILE__, util::WAR, ss.str());
                    throw std::runtime_error(ss.str());
                }
            }
            return result;
        }

        void validateRotamerTypes(const DihedralAngleDataTable& table, const ResidueLinkage& linkage)
        {
            for (auto& metadataVector : linkage.dihedralMetadata)
            {
                for (auto& metadata : metadataVector)
                {
                    if (table.entries[metadata].rotamer_type_ != linkage.rotamerType)
                    {
                        throw std::runtime_error(
                            "mismatching rotamer types in residue linkage: " + print(table, linkage));
                    }
                }
            }
        }

        void validateConformerMetadata(const DihedralAngleDataTable& table, const ResidueLinkage& linkage)
        {
            auto& dihedrals = linkage.rotatableBonds;
            auto& metadata = linkage.dihedralMetadata;
            size_t conformerCount = metadata[0].size();
            for (size_t n = 1; n < dihedrals.size(); n++)
            {
                if (metadata[n].size() != conformerCount)
                {
                    throw std::runtime_error(
                        "error: different number of conformers in linkage: " + print(table, linkage));
                }
            }

            for (size_t conformer = 0; conformer < conformerCount; conformer++)
            {
                double weight = table.entries[metadata[0][conformer]].weight_;
                for (size_t n = 1; n < dihedrals.size(); n++)
                {
                    if (std::fabs(weight - table.entries[metadata[n][conformer]].weight_) >= 1e-10)
                    {
                        throw std::runtime_error(
                            "error: conformers have different weight in linkage: " + print(table, linkage));
                    }
                }
            }
        }

        void throwMissingMetadataError(
            const ResidueLinkAttributes& names,
            const std::vector<DihedralAtoms>& dihedralAtoms,
            const std::vector<std::vector<size_t>>& metadataIndices)
        {
            size_t dihedralCount = dihedralAtoms.size();
            size_t missingCount = dihedralCount - metadataIndices.size();
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
                dihedralInfo.push_back(util::join(", ", atomNames));
            }
            stream << ": [" << util::join("], [", dihedralInfo) << "]";
            throw std::runtime_error(stream.str());
        }
    } // namespace

    unsigned long long generateResidueLinkageIndex()
    { // static keyword means it is created only once and persists beyond scope of code block.
        static unsigned long long s_ResidueLinkageIndex = 0;
        return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the
                                        // value in the copy
    }

    ResidueLink findResidueLink(std::pair<Residue*, Residue*> residues)
    {
        std::vector<Atom*> atoms;
        bool found = false;
        FindAtomsConnectingResidues(residues.first->getAtoms().at(0), residues.first, residues.second, &atoms, &found);
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

    void determineAtomsThatMove(std::vector<RotatableBond>& bonds)
    {
        for (auto& bond : bonds)
        {
            auto& atoms = bond.dihedralAtoms;
            std::vector<Atom*> atoms_that_move;
            atoms_that_move.push_back(atoms[2]);
            FindConnectedAtoms(atoms_that_move, atoms[1]);
            atoms_that_move.erase(atoms_that_move.begin());
            bond.movingAtoms = atoms_that_move;
        }
    }

    ResidueLinkage createResidueLinkage(const DihedralAngleDataTable& metadataTable, ResidueLink& link)
    {
        std::string firstId = residueStringId(link.residues.first);
        std::string secondId = residueStringId(link.residues.second);
        ResidueLinkAttributes linkAttributes = toAttributes(link);
        int local_debug = 1;
        if (local_debug > 0)
        {
            util::log(__LINE__, __FILE__, util::INF, "Maybe Finding connection between " + firstId + " :: " + secondId);
        }
        if (local_debug > 0)
        {
            util::log(__LINE__, __FILE__, util::INF, "Finding connection between " + firstId + " :: " + secondId);
            util::log(
                __LINE__,
                __FILE__,
                util::INF,
                "Connection atoms are from: " + link.atoms.first->getId() + " to " + link.atoms.second->getId());
        }
        std::vector<DihedralAtoms> dihedralAtoms = findRotatableDihedralsConnectingResidues(link);
        if (local_debug > 0)
        {
            util::log(__LINE__, __FILE__, util::INF, "Finding metadata for " + firstId + " :: " + secondId);
        }
        std::vector<std::vector<size_t>> metadata = findResidueLinkageMetadata(linkAttributes);
        if (local_debug > 0)
        {
            util::log(__LINE__, __FILE__, util::INF, "Metadata found:");
            for (auto& entry : metadata)
            {
                for (auto& dihedralAngleData : entry)
                {
                    const DihedralAngleData& entry = metadataTable.entries[dihedralAngleData];
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::INF,
                        util::join("_", {entry.atom1_, entry.atom2_, entry.atom3_, entry.atom4_}) + " : " +
                            util::join(
                                " ",
                                {std::to_string(entry.index_),
                                 entry.rotamer_name_,
                                 entry.dihedral_angle_name_,
                                 std::to_string(entry.default_angle)}));
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
            util::log(__LINE__, __FILE__, util::WAR, message);
        }
        // ensure that metadata and dihedral atoms have the same dimensions
        metadata.erase(metadata.begin() + dihedralAtoms.size(), metadata.end());

        auto& residues = link.residues;
        createHydrogenForPsiAngles(metadataTable, residues.second, dihedralAtoms, metadata);
        std::string linkageName = determineLinkageNameFromResidueNames(linkAttributes);
        std::vector<RotatableBond> bonds = createRotatableBonds(linkageName, dihedralAtoms, metadata);
        determineAtomsThatMove(bonds);

        unsigned long long index = generateResidueLinkageIndex();

        if (bonds.empty() || metadata[0].empty())
        {
            throw std::runtime_error("missing dihedrals or metadata in residue linkage: " + print(link));
        }

        const DihedralAngleData& firstMetadata = metadataTable.entries[metadata[0][0]];
        RotamerType rotamerType = firstMetadata.rotamer_type_;
        const std::vector<std::string>& cond1 = firstMetadata.residue1_conditions_;
        bool isDerivative = cond1.size() > 0 && cond1[0] == "derivative";
        ResidueLinkage linkage {link, bonds, metadata, rotamerType, isDerivative, index, linkageName};

        validateRotamerTypes(metadataTable, linkage);
        if (rotamerType == RotamerType::conformer)
        {
            validateConformerMetadata(metadataTable, linkage);
        }

        return linkage;
    }
} // namespace gmml
