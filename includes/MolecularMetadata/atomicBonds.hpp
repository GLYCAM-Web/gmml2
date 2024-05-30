#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include <string>
#include <utility>

namespace atomicBonds
{
    const double maxCutOff = 1.65;
    const double minCutOff = 0.7;
    // FUNCTIONS
    std::pair<double, double> getBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element);

    double getMaxBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element);

    template<class atomT> inline bool bondAtomsIfClose(atomT* atom1, atomT* atom2)
    {
        double maxLength = atomicBonds::getMaxBondLengthByAtomType(atom1->getElement(), atom2->getElement());
        if (withinDistance(maxLength, *atom1->getCoordinate(), *atom2->getCoordinate()))
        {
            atom1->addBond(atom2);
            return true;
        }
        return false;
    }

} // namespace atomicBonds

#endif /* INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP_ */
