#include "includes/CentralDataStructure/Selections/cyclePoints.hpp"

#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <algorithm>
#include <cstddef>

/* Below has following flaws:
  1) fail to find inner rotatable bonds here:
   __    __
__/  \__/  \__
  \__/  \__/
  If the above was within one residue.

  2) fail to find both connections of Res4 here:
 Res1 __
        \__ Res4
 Res2 __/

  3) It's ok with fused rings!

*/
std::vector<cds::Atom*> cdsSelections::FindCyclePoints(cds::Atom* atom, cds::Residue* residue)
{
    //    std::cout << "Entered FindCyclePoints with " << atom->getName() << std::endl;
    std::vector<cds::Atom*> rotation_points;
    std::vector<cds::Atom*> atom_path;
    bool found = false;
    if (residue->GetType() == cds::ResidueType::Protein)
    {
        cds::Atom* caAtom = residue->FindAtom("CA");
        // Find any cycle points. Note: starting from connection point to other residue
        cds::Atom* cycle_point;
        found = false;
        atom_path.clear();
        // This should only find cycles in Tyr, Trp etc. Not Asn, Ser, Thr as there aren't any unless bonding is messed
        // up.
        if (FindCyclePoint(atom, residue, atom, &atom_path, &found, cycle_point))
        {
            rotation_points.push_back(cycle_point);
            found = false;
            atom_path.clear();
            ClearAtomLabels(residue);
            // std::cout << "       >............." << std::endl;
            FindCyclePoint(caAtom, residue, caAtom, &atom_path, &found, cycle_point);
            rotation_points.push_back(cycle_point);
        }
        // Always want this at the end of the vector
        rotation_points.push_back(caAtom);
    }
    else
    {
        cds::Atom* rotation_point;
        atom_path.clear();
        found = false;
        // Find path to first cycle atom, i.e. anomeric carbon
        //       std::cout << "Non-protein, checking for cycles..." << std::endl;
        if (FindCyclePoint(atom, residue, atom, &atom_path, &found, rotation_point))
        {
            rotation_points.push_back(rotation_point);
        }
        // Ok, deal with non-protein non-cycles
        else
        { // Look for atom(s) with neighbors with no other neighbors within residue
            // I can't deal with non protein non cycles yet. How would I incode all the generic metadata?
            // Just set it all as rigid, and use the connecting atom as the "cycle point".
            rotation_points.push_back(atom);
        }
    }
    return rotation_points;
}

// Will not ignore fused rings. Explores everything to find all cycle points. Looks for cycle point closest to start
// atom.
bool cdsSelections::FindCyclePoint(cds::Atom* previous_atom, cds::Residue* residue, cds::Atom* current_atom,
                                   std::vector<cds::Atom*>* atom_path, bool* found_cycle_point, cds::Atom*& cycle_point)
{ // I definitely don't want cycles involving hydrogens (they only every form one bond, unless
    // bond by distance has bonded them).
    if (current_atom->getElement() != "H")
    {
        // Need this to explore everything. It will find same cycle point more than once, but that doesn't matter.
        current_atom->setLabels({"VisitedByFindCyclePoint"});
        atom_path->push_back(current_atom);
        std::vector<cds::Atom*> neighbors = current_atom->getNeighbors();
        for (std::vector<cds::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            cds::Atom* neighbor = *it1;
            // If not previous atom and not from a different residue
            if ((neighbor->getIndex() != previous_atom->getIndex()) && (residue->contains(neighbor)))
            {
                if (codeUtils::contains(*atom_path, neighbor)) // If we've been at this atom before
                {
                    if (*found_cycle_point) // If there are more than one cycle points
                    {
                        // Finds position of atoms in atom_path. Want earliest possible cycle point i.e. closest to
                        // start atom
                        std::ptrdiff_t new_cycle_position = std::distance(
                            atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), neighbor));
                        std::ptrdiff_t current_cycle_position = std::distance(
                            atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), cycle_point));
                        if (new_cycle_position < current_cycle_position)
                        {
                            cycle_point = neighbor;
                        }
                    }
                    else
                    {
                        *found_cycle_point = true;
                        cycle_point        = neighbor;
                    }
                }
                if (neighbor->getLabel() != "VisitedByFindCyclePoint") // Don't look back!
                {
                    FindCyclePoint(current_atom, residue, neighbor, atom_path, found_cycle_point, cycle_point);
                }
            }
        }
    }
    return *found_cycle_point;
}

// Find a neighbor of a cycle point to define the dihedral. Must not be in path used to come to this cycle_point.
// Oh gawd this code is horrible.
// The logic for selecting the atom to define a dihedral is messy.
// Comparing strings is messy when I really care about the number of e.g. C2 Vs O5, but need to factor in a C2 vs O
// comparision
cds::Atom* cdsSelections::FindCyclePointNeighbor(const std::vector<cds::Atom*> atom_path, cds::Atom* cycle_point,
                                                 cds::Residue* cyclePointResidue)
{
    cds::Atom* selected_neighbor;
    if (cycle_point->getName() == "CA") // This is a protein, and we always want the N atom.
    {
        selected_neighbor = cyclePointResidue->FindAtom("N");
    } // If this is a C2 like in Sia, then we always want the C1 atom unless that atom is in the linkage path (like
      // fructose 1-1)
    else if ((cycle_point->getName() == "C2") && (!codeUtils::contains(atom_path, cyclePointResidue->FindAtom("C1"))))
    {
        selected_neighbor = cyclePointResidue->FindAtom("C1");
    }
    else if (cycle_point->getName() == "C1")
    {
        selected_neighbor =
            cyclePointResidue->FindAtom("C2"); // We can always set phi to 180 this way, regardless of alpha/beta
    }
    else // This bit is overdone now, as I was looking for higher numbered atoms of C1, but now I know I always want C2,
         // so put that above.
    {
        std::vector<cds::Atom*> neighbors = cycle_point->getNeighbors();
        // Ok must first get a list of neighbors that weren't in the connection path
        std::vector<cds::Atom*> good_neighbors; // Couldn't think of a better name. Everybody needs these.
        for (auto& neighbor : neighbors)
        {
            if (!codeUtils::contains(atom_path, neighbor)) // If we've NOT been at this atom on way to cycle point
            {
                if (neighbor->getName().at(0) != 'H') // Don't find hydrogens. Later we swap out to use a hydrogen to
                                                      // define a dihedral, but that's a very specific one.
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        if (good_neighbors.size() == 0) // Ok take hydrogens then.
        {
            for (auto& neighbor : neighbors)
            {
                if (!codeUtils::contains(atom_path, neighbor)) // If we've NOT been at this atom on way to cycle point
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        selected_neighbor = good_neighbors.at(0); // Set to any to start. If there are not good_neighbors then you
                                                  // deserve to crash and burn std::cout << "Good neighbors are: ";
        for (std::vector<cds::Atom*>::iterator it1 = good_neighbors.begin(); it1 != good_neighbors.end(); ++it1)
        {
            cds::Atom* neighbor = *it1;
            if (selected_neighbor->getName().size() >= 2)
            {
                if (neighbor->getName().size() >=
                    2) // This is the only time I want to compare and select the larger number
                {
                    if (neighbor->getName().at(1) > selected_neighbor->getName().at(1))
                    {
                        selected_neighbor = neighbor;
                    }
                } // Otherwise any neighbor is ok. Yes I'm comparing char's, but that is fine unless C9 Vs C10, but in
                  // that case I don't care again.
            }
        }
    }
    return selected_neighbor;
}
