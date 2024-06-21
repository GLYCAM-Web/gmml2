#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp" //FindConnectedResidues()
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp" // For GetDescriptiveNameForGlycamResidueName in GetName function

#include <algorithm>
#include <memory>
#include <sstream>

using cds::ResidueLinkage;
using cds::RotatableDihedral;

namespace
{
    cds::ResidueLinkNames toNames(const cds::ResidueLink link)
    {
        return {
            {link.residues.first->getName(), link.residues.second->getName()},
            {   link.atoms.first->getName(),    link.atoms.second->getName()}
        };
    }

    std::string DetermineLinkageNameFromResidueNames(const cds::ResidueLinkNames link)
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

    //  This function splits that list into groups of 4 and creates RotatableDihedral object
    std::vector<cds::DihedralAtoms> SplitAtomVectorIntoRotatableDihedrals(bool isBranching,
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
                                        SplitAtomVectorIntoRotatableDihedrals(true, foundPath);
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
            SplitAtomVectorIntoRotatableDihedrals(false, connectingAtoms);
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

    std::vector<cds::RotatableDihedral> CreateRotatableDihedrals(const std::string& linkageName,
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

        std::vector<RotatableDihedral> rotatableDihedrals;
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

ResidueLinkage::ResidueLinkage(ResidueLink link) : link_(link)
{
    this->InitializeClass();
}

std::vector<RotatableDihedral> ResidueLinkage::GetRotatableDihedralsWithMultipleRotamers() const
{
    std::vector<RotatableDihedral> returningDihedrals;
    for (auto& entry : this->GetRotatableDihedrals())
    {
        if (entry.GetNumberOfRotamers() > 1)
        {
            returningDihedrals.push_back(entry);
        }
    }
    return returningDihedrals;
}

std::vector<RotatableDihedral>& ResidueLinkage::GetRotatableDihedralsRef()
{
    if (rotatableDihedrals_.empty())
    {
        std::stringstream ss;
        ss << "Error: RotatableDihedrals in this linkage is empty: " << link_.residues.first->getStringId() << "-"
           << link_.residues.second->getStringId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return rotatableDihedrals_;
}

std::vector<RotatableDihedral> ResidueLinkage::GetRotatableDihedrals() const
{
    if (rotatableDihedrals_.empty())
    {
        std::stringstream ss;
        ss << "Error: RotatableDihedrals in this linkage is empty: " << link_.residues.first->getStringId() << "-"
           << link_.residues.second->getStringId() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return rotatableDihedrals_;
}

int ResidueLinkage::GetNumberOfShapes(
    const bool likelyShapesOnly) const // Can have conformers (sets of rotamers) or permutations of rotamers
{
    int numberOfShapes = 1;
    if ((rotatableDihedrals_.empty()) || (rotatableDihedrals_.at(0).GetMetadata().empty()))
    {
        return numberOfShapes;
    }
    if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0)
    {
        for (auto& entry : rotatableDihedrals_)
        {
            numberOfShapes = (numberOfShapes * entry.GetNumberOfRotamers(likelyShapesOnly));
        }
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer") == 0)
    { // Conformer should mean that each dihedral will have the same number of metadata entries.
        // numberOfShapes = RotatableDihedrals_.size(); // This was correct for ASN for the wrong reason. 4 conformers
        // and 4 dihedrals...
        numberOfShapes = rotatableDihedrals_.at(0).GetNumberOfRotamers(likelyShapesOnly);
    }
    return numberOfShapes;
}

bool ResidueLinkage::CheckIfConformer() const
{
    if (rotatableDihedrals_.empty())
    {
        std::string errorMessage = "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.empty()\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().empty())
    {
        std::string errorMessage =
            "Error in ResidueLinkage::checkIfConformer as RotatableDihedrals_.at(0).GetMetadata().empty()\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    else
    {
        return (!(rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0));
    }
    return false; // Default to shut up the compiler. Shut up compiler gawd.
}

std::string ResidueLinkage::GetName() const
{
    if (!name_.empty())
    {
        return name_;
    }
    return DetermineLinkageNameFromResidueNames(toNames(link_));
}

void ResidueLinkage::AddNonReducingOverlapResidues(std::vector<cds::Residue*> extraResidues)
{
    nonReducingOverlapResidues_.insert(nonReducingOverlapResidues_.end(), extraResidues.begin(), extraResidues.end());
}

std::vector<cds::Residue*>& ResidueLinkage::GetNonReducingOverlapResidues()
{
    return nonReducingOverlapResidues_;
}

std::vector<cds::Residue*>& ResidueLinkage::GetReducingOverlapResidues()
{
    return reducingOverlapResidues_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void ResidueLinkage::SetDefaultShapeUsingMetadata()
{
    for (auto& entry : rotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(false, 0); // Default is first entry
    }
}

void ResidueLinkage::SetRandomShapeUsingMetadata()
{
    if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("permutation") == 0)
    {
        for (auto& entry : rotatableDihedrals_)
        {
            entry.SetRandomAngleEntryUsingMetadata();
        }
    }
    else if (rotatableDihedrals_.at(0).GetMetadata().at(0).rotamer_type_.compare("conformer") == 0)
    {
        int numberOfConformers = rotatableDihedrals_.at(0).GetMetadata().size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto& entry : rotatableDihedrals_)
        {
            entry.SetSpecificAngleEntryUsingMetadata(true, randomlySelectedConformerNumber);
        }
    }
}

void ResidueLinkage::SetSpecificShapeUsingMetadata(int shapeNumber)
{
    for (auto& entry : rotatableDihedrals_)
    {
        entry.SetSpecificAngleEntryUsingMetadata(false, shapeNumber);
    }
}

void ResidueLinkage::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
{
    for (auto& RotatableDihedral : rotatableDihedrals_)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (RotatableDihedral.SetSpecificShape(dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}

void ResidueLinkage::SetShapeToPrevious()
{
    for (typename std::vector<RotatableDihedral>::iterator RotatableDihedral = rotatableDihedrals_.begin();
         RotatableDihedral != rotatableDihedrals_.end(); ++RotatableDihedral)
    {
        RotatableDihedral->SetDihedralAngleToPrevious();
    }
    return;
}

void ResidueLinkage::DetermineAtomsThatMove()
{
    for (auto& dihedral : rotatableDihedrals_)
    {
        dihedral.DetermineAtomsThatMove();
    }
    return;
}

void ResidueLinkage::SimpleWiggleCurrentRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                                 std::vector<cds::Atom*>& overlapAtomSet2, const int angleIncrement)
{
    for (auto& RotatableDihedral : this->GetRotatableDihedrals())
    {
        RotatableDihedral.WiggleWithinCurrentRotamer(overlapAtomSet1, overlapAtomSet2, angleIncrement);
    }
}

void ResidueLinkage::SimpleWiggleCurrentRotamers(const std::array<std::vector<cds::Residue*>, 2>& residues,
                                                 const int angleIncrement)
{
    for (auto& dihedral : this->GetRotatableDihedrals())
    {
        const DihedralAngleDataVector rotamer {*dihedral.GetCurrentMetaData()};
        dihedral.WiggleUsingRotamers(rotamer, angleIncrement, residues);
    }
}

std::string ResidueLinkage::Print() const
{
    std::stringstream ss;
    ss << "ResidueLinkage Index: " << this->GetIndex() << ", Name: " << this->GetName()
       << ", NumberOfShapes: " << this->GetNumberOfShapes() << ", ids: " << this->GetFromThisResidue1()->getStringId()
       << "@" << link_.atoms.first->getName() << " -- " << this->GetToThisResidue2()->getStringId() << "@"
       << link_.atoms.second->getName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    for (auto& rotatableDihedral : this->GetRotatableDihedrals())
    {
        ss << rotatableDihedral.Print();
    }
    return ss.str();
}

// PRIVATE
void ResidueLinkage::InitializeClass()
{
    int local_debug = -1;
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Maybe Finding connection between " + link_.residues.first->getStringId() +
                      " :: " + link_.residues.second->getStringId());
    }
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Finding connection between " + link_.residues.first->getStringId() +
                      " :: " + link_.residues.second->getStringId());
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Connection atoms are from: " + link_.atoms.first->getId() + " to " + link_.atoms.second->getId());
    }
    std::vector<DihedralAtoms> dihedralAtoms = findRotatableDihedralsConnectingResidues(link_);
    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Finding metadata for " + link_.residues.first->getStringId() +
                      " :: " + link_.residues.second->getStringId());
    }
    DihedralAngleMetadata metadata = findResidueLinkageMetadata(toNames(link_));
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
    this->CreateHydrogenForPsiAngles(dihedralAtoms, metadata);
    rotatableDihedrals_ = CreateRotatableDihedrals(this->GetName(), dihedralAtoms, metadata);
    DetermineAtomsThatMove();
    this->SetIndex(this->GenerateIndex());
    this->DetermineResiduesForOverlapCheck(); // speedup overlap calcs
    return;
}

void ResidueLinkage::CreateHydrogenForPsiAngles(std::vector<DihedralAtoms>& dihedralAtoms,
                                                const DihedralAngleMetadata& metadata)
{
    for (size_t n = 0; n < dihedralAtoms.size(); n++)
    {
        for (auto& entry : metadata[n])
        {
            if (entry.dihedral_angle_name_ == "Psi" && entry.atom4_.at(0) == 'H')
            { // If it's a psi angle and is supposed to be defined by a H...
                Atom* atom     = dihedralAtoms[n].atoms[2];
                Atom* hydrogen = findHydrogenForPsiAngle(atom);
                if (hydrogen != nullptr)
                {
                    dihedralAtoms[n].atoms[3] = hydrogen;
                }
                else
                {
                    auto neighborsCoords = getCoordinatesFromAtoms(atom->getNeighbors());
                    Coordinate newCoord =
                        cds::CreateCoordinateForCenterAwayFromNeighbors(*atom->getCoordinate(), neighborsCoords);
                    std::unique_ptr<Atom> newAtom = std::make_unique<cds::Atom>("HHH", newCoord);
                    atom->addBond(newAtom.get());
                    link_.residues.second->addAtom(std::move(newAtom));
                }
            }
        }
    }
}

unsigned long long ResidueLinkage::GenerateIndex()
{ // static keyword means it is created only once and persists beyond scope of code block.
    static unsigned long long s_ResidueLinkageIndex = 0;
    return s_ResidueLinkageIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the
                                    // value in the copy
}

void ResidueLinkage::DetermineResiduesForOverlapCheck()
{
    reducingOverlapResidues_.clear();
    reducingOverlapResidues_.push_back(link_.residues.first);
    for (auto& neighbor : link_.residues.first->getNeighbors())
    {
        if (neighbor != link_.residues.second)
        {
            cdsSelections::FindConnectedResidues(reducingOverlapResidues_, neighbor);
        }
    }
    nonReducingOverlapResidues_.clear();
    nonReducingOverlapResidues_.push_back(link_.residues.second);
    for (auto& neighbor : link_.residues.second->getNeighbors())
    {
        if (neighbor != link_.residues.first)
        {
            cdsSelections::FindConnectedResidues(nonReducingOverlapResidues_, neighbor);
        }
    }
}
