#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/psiAngleHydrogen.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"

#include <sstream>

namespace
{
    cds::ResidueLinkNames toNames(const cds::ResidueLink link)
    {
        return {
            {link.residues.first->getName(), link.residues.second->getName()},
            {   link.atoms.first->getName(),    link.atoms.second->getName()}
        };
    }

    std::string determineLinkageNameFromResidueNames(const cds::ResidueLinkNames link)
    {
        std::string residue1Name = GlycamMetadata::GetDescriptiveNameForGlycamResidueName(link.residues.first);
        std::string residue2Name = GlycamMetadata::GetDescriptiveNameForGlycamResidueName(link.residues.second);
        std::string atom1Name    = link.atoms.first;
        std::string atom2Name    = link.atoms.second;
        char link1               = *atom1Name.rbegin(); //
        char link2               = *atom2Name.rbegin(); // Messy for Acetyl.
        std::stringstream linkageName;
        linkageName << residue1Name << link1 << "-" << link2 << residue2Name;
        return linkageName.str();
    }

    std::vector<cds::Residue*> connectedResidues(cds::Residue* residue, cds::Residue* block)
    {
        std::vector<cds::Residue*> result = {block};
        cdsSelections::FindConnectedResidues(result, residue);
        result.erase(result.begin());
        return result;
    }

    cds::DihedralAngleMetadata findResidueLinkageMetadata(cds::ResidueLinkNames link)
    {
        std::string firstAtom     = link.atoms.first;
        std::string secondAtom    = link.atoms.second;
        std::string firstResidue  = link.residues.first;
        std::string secondResidue = link.residues.second;
        cds::DihedralAngleMetadata matching_entries =
            gmml::MolecularMetadata::GLYCAM::getDihedralAngleDataEntriesForLinkage(firstAtom, firstResidue, secondAtom,
                                                                                   secondResidue);
        if (matching_entries.empty())
        {
            matching_entries = gmml::MolecularMetadata::GLYCAM::getDihedralAngleDataEntriesForLinkage(
                secondAtom, secondResidue, firstAtom, firstResidue);
        }
        if (matching_entries.empty())
        {
            std::stringstream ss;
            ss << "No Metadata entries found for connection between " << firstResidue << "@" << firstAtom << " and "
               << secondResidue << "@" << secondAtom << "\n";
            ss << "Note that order should be reducing atom - anomeric atom, but I've tried reversing the order and it "
                  "didn't fix the issue.\n";
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            throw std::runtime_error(ss.str());
        }
        return matching_entries;
    }

    std::vector<cds::RotatableDihedral> createRotatableDihedrals(const std::string& linkageName,
                                                                 const std::vector<cds::DihedralAtoms>& dihedralAtoms,
                                                                 const cds::DihedralAngleMetadata& metadata)
    {
        if (metadata.size() > dihedralAtoms.size())
        {
            std::string message =
                "Found metadata for rotatable bonds that do not exist.\nCheck both dihedralangledata metadata "
                "and ResidueLinkage::FindRotatableDihedralsConnectingResidues.\nNote this is normal for a sialic acid "
                "with multiple 2-7, 2-8 and or 2-9 linkages and this warning can be ignored\n";
            gmml::log(__LINE__, __FILE__, gmml::WAR, message);
        }

        std::vector<cds::RotatableDihedral> rotatableDihedrals;
        rotatableDihedrals.reserve(dihedralAtoms.size());
        for (size_t n = 0; n < dihedralAtoms.size(); n++)
        {
            auto& currentMetadata = metadata[n];
            if (!currentMetadata.empty())
            {
                rotatableDihedrals.emplace_back(
                    cds::RotatableDihedral {dihedralAtoms[n].isBranching, dihedralAtoms[n].atoms, currentMetadata});
            }
            else
            {
                std::stringstream ss;
                ss << "Problem with the metadata found in gmml for this linkage. No metadata found for dihedral with "
                      "bond: ";
                ss << linkageName << "\n";
                ss << "At index: " << n << "\n";
                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                throw std::runtime_error(ss.str());
            }
        }
        return rotatableDihedrals;
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
        dihedral.movingCoordinates = getCoordinatesFromAtoms(atoms_that_move);
    }
}

void cds::determineResiduesForOverlapCheck(ResidueLinkage& linkage)
{
    auto& residues                     = linkage.link.residues;
    linkage.reducingOverlapResidues    = connectedResidues(residues.first, residues.second);
    linkage.nonReducingOverlapResidues = connectedResidues(residues.second, residues.first);
}

cds::ResidueLinkage cds::createResidueLinkage(ResidueLink& link)
{
    int local_debug = -1;
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Maybe Finding connection between " + link.residues.first->getStringId() +
                      " :: " + link.residues.second->getStringId());
    }
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Finding connection between " + link.residues.first->getStringId() +
                      " :: " + link.residues.second->getStringId());
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Connection atoms are from: " + link.atoms.first->getId() + " to " + link.atoms.second->getId());
    }
    std::vector<DihedralAtoms> dihedralAtoms = cdsSelections::findRotatableDihedralsConnectingResidues(link);
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Finding metadata for " + link.residues.first->getStringId() +
                      " :: " + link.residues.second->getStringId());
    }
    DihedralAngleMetadata metadata = findResidueLinkageMetadata(toNames(link));
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Metadata found:");
        for (auto& entry : metadata)
        {
            for (auto& dihedralAngleData : entry)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, dihedralAngleData.print());
            }
        }
    }
    auto& residues = link.residues;
    createHydrogenForPsiAngles(residues.second, dihedralAtoms, metadata);
    std::string name                         = determineLinkageNameFromResidueNames(toNames(link));
    std::vector<RotatableDihedral> dihedrals = createRotatableDihedrals(name, dihedralAtoms, metadata);
    determineAtomsThatMove(dihedrals);

    unsigned long long index        = generateResidueLinkageIndex();
    auto reducingOverlapResidues    = connectedResidues(residues.first, residues.second);
    auto nonReducingOverlapResidues = connectedResidues(residues.second, residues.first);

    if (dihedrals.empty() || dihedrals[0].metadataVector.empty())
    {
        throw std::runtime_error("missing dihedrals or metadata in residue linkage");
    }

    auto rotamerType = dihedrals[0].metadataVector[0].rotamer_type_;
    for (auto& dihedral : dihedrals)
    {
        for (auto& metadata : dihedral.metadataVector)
        {
            if (metadata.rotamer_type_ != rotamerType)
            {
                throw std::runtime_error("mismatching rotamer types in residue linkage");
            }
        }
    }

    return ResidueLinkage(link, dihedrals, rotamerType, index, name, reducingOverlapResidues,
                          nonReducingOverlapResidues);
}
