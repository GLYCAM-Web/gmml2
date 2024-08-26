#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICBONDING_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICBONDING_HPP

#include "includes/CentralDataStructure/atom.hpp"

namespace cds
{
    bool isWithinBondingDistance(const Atom* atom, const Atom* otherAtom);
    void addBond(Atom* atom, Atom* otherAtom);
    bool bondAtomsIfClose(cds::Atom* atom1, cds::Atom* atom2);
} // namespace cds

#endif
