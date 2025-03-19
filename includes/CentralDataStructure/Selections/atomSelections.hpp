#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include <vector>
#include <string>

using cds::Atom;
using cds::Residue;

namespace cdsSelections
{
    Atom* getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string& queryNumber);
    void FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom);
    Atom* selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms);
    std::vector<Atom*> findCycleAtoms(cds::Atom* const starterAtom);
    Atom* guessAnomericAtomByInternalNeighbors(const std::vector<cds::Atom*> atoms);
    std::vector<Atom*> FindHeavyAtoms(std::vector<Atom*> queryAtoms);
    unsigned long int CountAtomsWithinBondingDistance(const Atom* queryAtom, std::vector<Atom*> otherAtoms);
} // namespace cdsSelections
#endif
