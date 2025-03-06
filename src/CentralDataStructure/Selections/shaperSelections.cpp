#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"

#include "includes/CentralDataStructure/Selections/cyclePoints.hpp"
#include "includes/CentralDataStructure/Selections/linkageBranches.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <sstream>
#include <vector>

bool cdsSelections::FindPathBetweenTwoAtoms(cds::Atom* current_atom, cds::Residue* currentResidue,
                                            cds::Atom* target_atom, cds::Residue* targetResidue,
                                            std::vector<cds::Atom*>* atom_path, bool* found)
{
    current_atom->setLabels({"VistedByFindPathBetweenTwoAtoms"});
    std::vector<cds::Atom*> neighbors = current_atom->getNeighbors();
    for (std::vector<cds::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        cds::Atom* neighbor = *it1;
        if (neighbor->getIndex() == target_atom->getIndex())
        {
            *found = true;
            atom_path->push_back(neighbor);
        }
        // If not found && not previously visited atom && ( if neighbor residue is current residue || target_atom
        // residue)
        if ((*found == false) && (neighbor->getLabel() != "VistedByFindPathBetweenTwoAtoms") &&
            ((currentResidue->contains(neighbor)) || (targetResidue->contains(neighbor))))
        {
            FindPathBetweenTwoAtoms(neighbor, currentResidue, target_atom, targetResidue, atom_path, found);
        }
    }
    if (*found) // As you fall back out from the found target, create a list of the atoms.
    {
        atom_path->push_back(current_atom);
    }
    return *found;
}

// I want a generic recursive function, where I can pass in the condition(s). Lots of Repeating code here.
// This one was written before the others. Could update with previous atom being passed in, though that makes the
// initial call confusing...
void cdsSelections::FindAtomsConnectingResidues(cds::Atom* current_atom, const cds::Residue* currentResidue,
                                                const cds::Residue* otherResidue,
                                                std::vector<cds::Atom*>* connecting_atoms, bool* found_neighbor)
{
    current_atom->setLabels({"VisitedByFindAtomsConnectingResidues"});
    for (auto& neighbor : current_atom->getNeighbors())
    {
        if (otherResidue->contains(neighbor))
        {
            *found_neighbor = true;
            connecting_atoms->push_back(current_atom);
            connecting_atoms->push_back(neighbor);
        }
        // If haven't visited this atom already AND don't move onto other residues
        else if ((neighbor->getLabel() != "VisitedByFindAtomsConnectingResidues") &&
                 (currentResidue->contains(neighbor)))
        {
            FindAtomsConnectingResidues(neighbor, currentResidue, otherResidue, connecting_atoms, found_neighbor);
        }
    }
    return;
}

void cdsSelections::ClearAtomLabels(cds::Residue* residue)
{
    for (auto& atom : residue->getAtoms())
    {
        atom->clearLabels();
    }
    return;
}

cds::ResidueLinkage* cdsSelections::selectLinkageWithIndex(std::vector<cds::ResidueLinkage>& inputLinkages,
                                                           const long long unsigned int indexQuery)
{
    for (auto& linkage : inputLinkages)
    {
        if (linkage.index == indexQuery)
        {
            return &linkage;
        }
    }
    // Error
    std::stringstream ss;
    ss << "Linkage numbered " << indexQuery << " not found in linkages for this carbohydrate\n";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    throw std::runtime_error(ss.str());
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
/* This is a straight copy from glycoprotein_builder. I need a high level class that deals with both
 * cds::ResidueLinkages, ring shapes etc. That way I can create X shapes of a molecule. For now this will do to figure
 * out some implementation details like file naming.
 */
std::vector<cds::ResidueLinkage>
cdsSelections::SplitLinkagesIntoPermutants(std::vector<cds::ResidueLinkage>& inputLinkages)
{
    std::vector<cds::ResidueLinkage> sortedLinkages;
    for (auto& linkage : inputLinkages)
    {
        if (linkage.rotamerType == GlycamMetadata::RotamerType::conformer)
        {
            sortedLinkages.push_back(linkage);
        }
        else // if not a conformer
        {
            std::vector<size_t> rotatableDihedralIndices = cds::rotatableDihedralsWithMultipleRotamers(
                linkage.dihedralMetadata); // only want the rotatabe dihedrals within a linkage
                                           // that have multiple rotamers. Some bonds won't.
            for (size_t index : rotatableDihedralIndices)
            {
                cds::ResidueLinkage splitLinkage = linkage; // Copy it to get correct info into class
                splitLinkage.rotatableDihedrals  = {linkage.rotatableDihedrals[index]};
                splitLinkage.dihedralMetadata    = {linkage.dihedralMetadata[index]};
                sortedLinkages.push_back(splitLinkage);
            }
        }
    }
    return sortedLinkages;
}

//  This function splits that list into groups of 4 and creates RotatableDihedral object
std::vector<cds::DihedralAtoms>
cdsSelections::splitAtomVectorIntoRotatableDihedrals(const std::vector<cds::Atom*>& atoms)
{
    // Ok looking for sets of four atoms, but shifting along vector by one atom for each dihedral.
    //  So four atoms will make one rotatable bond, five will make two bonds, six will make three etc.
    if (atoms.size() < 4)
    {
        std::stringstream ss;
        ss << "ERROR in splitAtomVectorIntoRotatableDihedrals, not enough atoms in atom vector: " << atoms.size()
           << "\n";
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
            dihedrals.emplace_back(cds::DihedralAtoms {atoms[n], atoms[n + 1], atoms[n + 2], atoms[n + 3]});
        }
        return dihedrals;
    }
}

//  generates a list of linearly connected atoms that define the rotatable bonds
std::vector<cds::DihedralAtoms> cdsSelections::findRotatableDihedralsConnectingResidues(const cds::ResidueLink& link)
{
    // Going to ignore tags etc.
    // Given two residues that are connected. Find connecting atoms.
    // Search neighbors other than connected atom. Ie search out in both directions, but remain within same residue.
    // Warning, residue may have fused cycles!
    // Will fail for non-protein residues without cycles. As don't have a non-rotatable bond to anchor from. Can
    // code that later (and deal with branches from these residues).

    std::vector<cds::Atom*> firstResidueCyclePoints  = FindCyclePoints(link.atoms.first, link.residues.first);
    //    std::cout << "Moving onto second residue.\n";
    std::vector<cds::Atom*> secondResidueCyclePoints = FindCyclePoints(link.atoms.second, link.residues.second);
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
    std::vector<cds::DihedralAtoms> RotatableDihedrals = splitAtomVectorIntoRotatableDihedrals(connectingAtoms);
    // Add any linkage branches (in 2-7 and 2-8) to the rest.
    RotatableDihedrals.insert(RotatableDihedrals.end(), rotatableDihedralsInBranches.begin(),
                              rotatableDihedralsInBranches.end());
    return RotatableDihedrals;
}
