#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include <string>
#include <utility>

namespace atomicBonds
{
    enum AtomicElement
    {
        C       = 6,
        N       = 7,
        O       = 8,
        P       = 15,
        S       = 16,
        Unknown = 0
    };

    const double maxCutOff = 1.65;
    const double minCutOff = 0.7;

    AtomicElement toElement(const std::string& str);
    std::pair<double, double> getBondLengthByAtomType(AtomicElement atom1Element, AtomicElement atom2Element);

    double getMaxBondLengthByAtomType(AtomicElement atom1Element, AtomicElement atom2Element);

    template<class atomT> inline bool bondAtomsIfClose(atomT* atom1, atomT* atom2)
    {
        double maxLength =
            atomicBonds::getMaxBondLengthByAtomType(toElement(atom1->getElement()), toElement(atom2->getElement()));
        if (withinDistance(maxLength, *atom1->getCoordinate(), *atom2->getCoordinate()))
        {
            atom1->addBond(atom2);
            return true;
        }
        return false;
    }

} // namespace atomicBonds

#endif
