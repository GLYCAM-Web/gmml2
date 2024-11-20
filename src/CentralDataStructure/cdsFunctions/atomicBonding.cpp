#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/MolecularMetadata/elements.hpp"

bool cds::isWithinBondingDistance(const Atom* atom, const Atom* otherAtom)
{
    double maxLength = MolecularMetadata::maxBondLengthByAtomType(
        MolecularMetadata::toElement(atom->getElement()), MolecularMetadata::toElement(otherAtom->getElement()));
    return withinDistance(maxLength, atom->coordinate(), otherAtom->coordinate());
}

void cds::addBond(Atom* atom, Atom* otherAtom)
{
    atom->addNeighbor("atomicBond", otherAtom);
}
