#include "includes/CentralDataStructure/Selections/linkageBranches.hpp"

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/Selections/cyclePoints.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <tuple>
#include <sstream>

// After getting connection atoms between two residues, we have the linear path between them e.g.:
// (Residue1 C2 - O7 - C7 - O6 - C6 Residue2), where C2 and C6 are the cycle points. In cases like 2-7 or 2-8 linkages,
// we have significant branches from this linear path
//                    /
//                   C8-O8
//                  /
//                 C9-O9
// And we need to set reasonable values for the C9-C8, and C8-C7 dihedral angles.
// This is further complicated by the possibility of "DNeup5Aca2-7[DNeup5Aca2-8]DNeup5Aca2-OH", where one of the
// branches is part of another linkage. This will work for 7/8/9 linked sialic acids, but if we get more heavily
// branched linkages this will need to change to iteratively find branches from branches
//

void cdsSelections::FindEndsOfBranchesFromLinkageAtom(cds::Atom* currentAtom, cds::Atom* previousAtom,
                                                      cds::Residue* residue, Branch* branch)
{
    branch->ChangeDepth(1);
    currentAtom->setLabels({"VistedByFindEndsOfBranchesFromLinkageAtom"});
    bool deadEndAtom              = true;
    bool connectsToAnotherResidue = false;
    for (auto& neighbor : currentAtom->getNeighbors())
    {
        if (neighbor->getLabel() != "VistedByFindEndsOfBranchesFromLinkageAtom" && *neighbor != *previousAtom &&
            residue->contains(neighbor)) // Don't explore across residues.
        {
            if (neighbor->getNeighbors().size() > 1)
            {
                deadEndAtom = false;
                FindEndsOfBranchesFromLinkageAtom(neighbor, currentAtom, residue, branch);
                branch->ChangeDepth(-1);
            }
        }
        if (!residue->contains(neighbor))
        {
            connectsToAnotherResidue = true;
        }
    }
    if (deadEndAtom && !connectsToAnotherResidue && branch->GetDepth() > 1 && branch->AtMaxDepth())
    {
        branch->SetEnd(currentAtom);
    }
    return;
}

std::tuple<std::vector<cds::DihedralAtoms>, std::vector<cds::Atom*>>
cdsSelections::findRotatableDihedralsinBranchesConnectingResidues(const cds::ResidueLink& link,
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
        cds::Atom* neighbor1 = cdsSelections::FindCyclePointNeighbor(connectingAtoms, cyclePoint1, link.residues.first);
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
                        if (!codeUtils::contains(connectingAtoms, neighbor)) // if not in the vector
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
                                cds::Atom* neighbor = cdsSelections::FindCyclePointNeighbor(foundPath, branch.GetRoot(),
                                                                                            link.residues.second);
                                foundPath.push_back(neighbor);
                                std::vector<cds::DihedralAtoms> temp =
                                    splitAtomVectorIntoRotatableDihedrals(true, foundPath);
                                codeUtils::insertInto(rotatableDihedralsInBranches, temp);
                            }
                        }
                    }
                }
            }
        } // End dealing with branching linkages
    }
    return {rotatableDihedralsInBranches, connectingAtoms};
}
