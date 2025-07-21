#ifndef INCLUDES_SELECTIONS_CYCLEPOINTS_HPP
#define INCLUDES_SELECTIONS_CYCLEPOINTS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"

#include <vector>

namespace gmml
{
    std::vector<Atom*> FindCyclePoints(Atom* atom, Residue* residue);

    bool FindCyclePoint(
        Atom* previous_atom,
        Residue* residue,
        Atom* current_atom,
        std::vector<Atom*>* atom_path,
        bool* found_cycle_point,
        Atom*& cycle_point);

    Atom* FindCyclePointNeighbor(const std::vector<Atom*> atom_path, Atom* cycle_point, Residue* cyclePointResidue);
} // namespace gmml

#endif
