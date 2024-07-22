#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_CYCLEPOINTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_CYCLEPOINTS_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include <vector>

namespace cdsSelections
{
    std::vector<cds::Atom*> FindCyclePoints(cds::Atom* atom, cds::Residue* residue);
    bool FindCyclePoint(cds::Atom* previous_atom, cds::Residue* residue, cds::Atom* current_atom,
                        std::vector<cds::Atom*>* atom_path, bool* found_cycle_point, cds::Atom*& cycle_point);
    cds::Atom* FindCyclePointNeighbor(const std::vector<cds::Atom*> atom_path, cds::Atom* cycle_point,
                                      cds::Residue* cyclePointResidue);
} // namespace cdsSelections
#endif
