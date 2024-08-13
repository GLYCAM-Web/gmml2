#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ATOMTOCOORDINATEINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_ATOMTOCOORDINATEINTERFACE_HPP

// This file replaces some of the functions that are in geometrytopology, the ones that take MoleculeModeling classes,
// and replaces them It still calls the Coordinate accepting classes in geometrytopology. The idea is to separate
// geometrytopology into the parts that use coordinates only.
#include "includes/CentralDataStructure/atom.hpp"

using cds::Coordinate;

namespace cds
{
    void moveConnectedAtomsAccordingToBondLength(cds::Atom* a, cds::Atom* b);
} // namespace cds
#endif
