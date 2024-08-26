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

} // namespace MolecularMetadata

#endif
