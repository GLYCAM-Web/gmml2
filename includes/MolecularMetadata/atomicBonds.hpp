#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include "includes/MolecularMetadata/elements.hpp"

#include <string>
#include <utility>

namespace MolecularMetadata
{
    const double maxCutOff = 1.65;
    const double minCutOff = 0.7;

    std::pair<double, double> getBondLengthByAtomType(Element atom1Element, Element atom2Element);

    double getMaxBondLengthByAtomType(Element atom1Element, Element atom2Element);

    template<class atomT> inline bool bondAtomsIfClose(atomT* atom1, atomT* atom2)
    {
        double maxLength = MolecularMetadata::getMaxBondLengthByAtomType(toElement(atom1->getElement()),
                                                                         toElement(atom2->getElement()));
        if (withinDistance(maxLength, *atom1->getCoordinate(), *atom2->getCoordinate()))
        {
            atom1->addBond(atom2);
            return true;
        }
        return false;
    }

} // namespace MolecularMetadata

#endif
