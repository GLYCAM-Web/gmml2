#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"

#include <algorithm>
#include <functional>
#include <memory>
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

    //  This function splits that list into groups of 4 and creates RotatableDihedral object
    std::vector<cds::DihedralAtoms> splitAtomVectorIntoRotatableDihedrals(bool isBranching,
                                                                          const std::vector<cds::Atom*>& atoms)
    {
        // Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
        //  So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
        if (atoms.size() < 4)
        {
            std::stringstream ss;
            ss << "ERROR in ResidueLinkage::SplitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: "
               << atoms.size() << "\n";
            ss << "This should be 4 or something is very wrong\n";
            ss << "If there are atoms, here are the ids:\n";
            for (auto& atom : atoms)
            {
                ss << atom->getId() << "\n";
            }
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            throw std::runtime_error(ss.str());
        }
        else
        {
            std::vector<cds::DihedralAtoms> dihedrals;
            for (size_t n = 0; n < atoms.size() - 3; n++)
            {
                dihedrals.emplace_back(cds::DihedralAtoms {
                    isBranching, {atoms[n], atoms[n + 1], atoms[n + 2], atoms[n + 3]}
                });
            }
            return dihedrals;
        }
    }

    std::tuple<std::vector<cds::DihedralAtoms>, std::vector<cds::Atom*>>
    findRotatableDihedralsinBranchesConnectingResidues(const cds::ResidueLink& link,
                                                       const std::vector<cds::Atom*>& residueCyclePoints)
    {
        std::vector<cds::DihedralAtoms> rotatableDihedralsInBranches;
        std::vector<cds::Atom*> connectingAtoms;
        for (long unsigned int i = 0; i < residueCyclePoints.size(); i = i + 2)
        {
            cds::Atom* cyclePoint1 = residueCyclePoints.at(i);
            cds::Atom* cyclePoint2 = residueCyclePoints.at(i + 1);

            bool found = false;
            connectingAtoms.clear();
            cdsSelections::FindPathBetweenTwoAtoms(cyclePoint1, link.residues.first, cyclePoint2, link.residues.second,
                                                   &connectingAtoms, &found);
            cdsSelections::ClearAtomLabels(
                link.residues.first); // ToDo change to free function or member function that clears labels.
            cdsSelections::ClearAtomLabels(link.residues.second);
            // Find neighboring atoms needed to define dihedral. Pass in connecting atoms so don't find any of those.
            cds::Atom* neighbor1 =
                cdsSelections::FindCyclePointNeighbor(connectingAtoms, cyclePoint1, link.residues.first);
            cds::Atom* neighbor2 =
                cdsSelections::FindCyclePointNeighbor(connectingAtoms, cyclePoint2, link.residues.second);
            // Insert these neighbors into list of connecting atoms, at beginning and end of vector.
            // connecting_atoms gets populated as it falls out, so list is reversed from what you'd expect
            std::reverse(connectingAtoms.begin(), connectingAtoms.end());
            connectingAtoms.insert(connectingAtoms.begin(), neighbor1);
            connectingAtoms.push_back(neighbor2);
            cdsSelections::ClearAtomLabels(
                link.residues.first); // ToDo change to free function or member function that clears labels.
            cdsSelections::ClearAtomLabels(link.residues.second);
            // This mess was made to address the branching in 2-7 and 2-8 linkages.
            // These branches are long enough that they need default torsions set.
            if (connectingAtoms.size() > 4) // Otherwise there are no torsions
            {                               // Only viable linkages. Throw if not >4?
                for (typename std::vector<cds::Atom*>::iterator it = connectingAtoms.begin() + 1;
                     it != connectingAtoms.end() - 1; ++it)
                { //
                    cds::Atom* connectionAtom = (*it);
                    if ((connectionAtom != cyclePoint1) && (connectionAtom != cyclePoint2))
                    {
                        for (auto& neighbor : connectionAtom->getNeighbors()) // Need an interator
                        {
                            if (std::find(connectingAtoms.begin(), connectingAtoms.end(), neighbor) ==
                                connectingAtoms.end()) // if not in the vector
                            {
                                cdsSelections::Branch branch(connectionAtom);
                                cdsSelections::FindEndsOfBranchesFromLinkageAtom(neighbor, connectionAtom,
                                                                                 link.residues.second, &branch);
                                if (branch.IsBranchFound())
                                {
                                    found = false;
                                    std::vector<cds::Atom*> foundPath;
                                    // This fills in foundPath:
                                    cdsSelections::FindPathBetweenTwoAtoms(branch.GetRoot(), link.residues.first,
                                                                           branch.GetEnd(), link.residues.second,
                                                                           &foundPath, &found);
                                    cds::Atom* neighbor = cdsSelections::FindCyclePointNeighbor(
                                        foundPath, branch.GetRoot(), link.residues.second);
                                    foundPath.push_back(neighbor);
                                    std::vector<cds::DihedralAtoms> temp =
                                        splitAtomVectorIntoRotatableDihedrals(true, foundPath);
                                    rotatableDihedralsInBranches.insert(rotatableDihedralsInBranches.end(),
                                                                        temp.begin(), temp.end());
                                }
                            }
                        }
                    }
                }
            } // End dealing with branching linkages
        }
        return {rotatableDihedralsInBranches, connectingAtoms};
    }

    //  generates a list of linearly connected atoms that define the rotatable bonds
    std::vector<cds::DihedralAtoms> findRotatableDihedralsConnectingResidues(const cds::ResidueLink& link)
    {
        // Going to ignore tags etc.
        // Given two residues that are connected. Find connecting atoms.
        // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
        // Warning, residue may have fused cycles!
        // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can
        // code that later (and deal with branches from these residues).

        std::vector<cds::Atom*> firstResidueCyclePoints =
            cdsSelections::FindCyclePoints(link.atoms.first, link.residues.first);
        //    std::cout << "Moving onto second residue.\n";
        std::vector<cds::Atom*> secondResidueCyclePoints =
            cdsSelections::FindCyclePoints(link.atoms.second, link.residues.second);
        // Need to reverse one of these, so when concatenated, they are ordered ok. This might not be ok.
        // std::reverse(to_this_residue2_cycle_points.begin(), to_this_residue2_cycle_points.end());
        std::reverse(firstResidueCyclePoints.begin(), firstResidueCyclePoints.end());
        // Now concatenate:
        firstResidueCyclePoints.insert(firstResidueCyclePoints.end(), secondResidueCyclePoints.begin(),
                                       secondResidueCyclePoints.end());
        // Now that have a list of rotation points. Split into pairs and find rotatable bonds between them
        auto branchResult = findRotatableDihedralsinBranchesConnectingResidues(link, firstResidueCyclePoints);
        std::vector<cds::DihedralAtoms>& rotatableDihedralsInBranches = std::get<0>(branchResult);
        std::vector<cds::Atom*>& connectingAtoms                      = std::get<1>(branchResult);
        std::vector<cds::DihedralAtoms> RotatableDihedrals =
            splitAtomVectorIntoRotatableDihedrals(false, connectingAtoms);
        // Add any linkage branches (in 2-7 and 2-8) to the rest.
        RotatableDihedrals.insert(RotatableDihedrals.end(), rotatableDihedralsInBranches.begin(),
                                  rotatableDihedralsInBranches.end());
        return RotatableDihedrals;
    }

    cds::DihedralAngleMetadata findResidueLinkageMetadata(cds::ResidueLinkNames link)
    {
        std::string firstAtom     = link.atoms.first;
        std::string secondAtom    = link.atoms.second;
        std::string firstResidue  = link.residues.first;
        std::string secondResidue = link.residues.second;
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer DihedralAngleMetadata;
        cds::DihedralAngleMetadata matching_entries =
            DihedralAngleMetadata.GetEntriesForLinkage(firstAtom, firstResidue, secondAtom, secondResidue);
        if (matching_entries.empty())
        {
            matching_entries =
                DihedralAngleMetadata.GetEntriesForLinkage(secondAtom, secondResidue, firstAtom, firstResidue);
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

    cds::Atom* findHydrogenForPsiAngle(const cds::Atom* atom)
    {
        for (auto& neighbor : atom->getNeighbors())
        {
            if (neighbor->getName().at(0) == 'H')
            {
                return neighbor;
            }
        }
        return nullptr;
    }

    void createHydrogenForPsiAngles(cds::Residue* residue, std::vector<cds::DihedralAtoms>& dihedralAtoms,
                                    const cds::DihedralAngleMetadata& metadata)
    {
        for (size_t n = 0; n < dihedralAtoms.size(); n++)
        {
            for (auto& entry : metadata[n])
            {
                if (entry.dihedral_angle_name_ == "Psi" && entry.atom4_.at(0) == 'H')
                { // If it's a psi angle and is supposed to be defined by a H...
                    cds::Atom* atom     = dihedralAtoms[n].atoms[2];
                    cds::Atom* hydrogen = findHydrogenForPsiAngle(atom);
                    if (hydrogen != nullptr)
                    {
                        dihedralAtoms[n].atoms[3] = hydrogen;
                    }
                    else
                    {
                        Coordinate newCoord = cds::coordinateOppositeToNeighborAverage(
                            *atom->getCoordinate(), cds::getCoordinatesFromAtoms(atom->getNeighbors()), 1.0);
                        std::unique_ptr<cds::Atom> newAtom = std::make_unique<cds::Atom>("HHH", newCoord);
                        atom->addBond(newAtom.get());
                        residue->addAtom(std::move(newAtom));
                    }
                }
            }
        }
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
    std::vector<DihedralAtoms> dihedralAtoms = findRotatableDihedralsConnectingResidues(link);
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

    return ResidueLinkage(link, dihedrals, index, name, reducingOverlapResidues, nonReducingOverlapResidues);
}