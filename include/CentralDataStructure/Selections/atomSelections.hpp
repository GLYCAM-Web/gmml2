#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_ATOMSELECTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/elements.hpp"

#include <string>
#include <vector>

namespace gmml
{
    Atom* getNonCarbonHeavyAtomNumbered(std::vector<Atom*> atoms, const std::string& queryNumber);
    void FindConnectedAtoms(std::vector<Atom*>& visitedAtoms, Atom* currentAtom);
    Atom* selectNeighborNotInAtomVector(const Atom* atomWithNeighbors, std::vector<Atom*> queryAtoms);
    std::vector<Atom*> findCycleAtoms(Atom* const starterAtom);
    Atom* guessAnomericAtomByInternalNeighbors(const std::vector<Atom*> atoms);
    std::vector<size_t> FindHeavyAtoms(const std::vector<Element>& elements);
} // namespace gmml

#endif
