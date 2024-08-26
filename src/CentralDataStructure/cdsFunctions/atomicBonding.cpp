#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/Geometry/functions.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/MolecularMetadata/elements.hpp"

bool cds::isWithinBondingDistance(const Atom* atom, const Atom* otherAtom)
{
    double maxLength = MolecularMetadata::getMaxBondLengthByAtomType(
        MolecularMetadata::toElement(atom->getElement()), MolecularMetadata::toElement(otherAtom->getElement()));
    return withinDistance(maxLength, *atom->getCoordinate(), *otherAtom->getCoordinate());
}

void cds::addBond(Atom* atom, Atom* otherAtom)
{
    atom->addNeighbor("atomicBond", otherAtom);
}

bool cds::bondAtomsIfClose(cds::Atom* atom1, cds::Atom* atom2)
{
    double maxLength = MolecularMetadata::getMaxBondLengthByAtomType(MolecularMetadata::toElement(atom1->getElement()),
                                                                     MolecularMetadata::toElement(atom2->getElement()));
    if (withinDistance(maxLength, *atom1->getCoordinate(), *atom2->getCoordinate()))
    {
        addBond(atom1, atom2);
        return true;
    }
    return false;
}
