#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include "includes/MolecularMetadata/elements.hpp"

#include <string>

namespace MolecularMetadata
{
    double hydrogenCovalentBondMaxLength();
    double maxBondLengthByAtomType(Element atom1Element, Element atom2Element);
    double specificBondLength(const std::string& query1, const std::string& query2);
    double deoxyHydrogenBondLength();

} // namespace MolecularMetadata

#endif
