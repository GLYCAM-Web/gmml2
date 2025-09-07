#ifndef INCLUDE_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    Atom* getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string& queryNumber);
    void FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom);
    Atom* selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms);
} // namespace gmml

#endif
