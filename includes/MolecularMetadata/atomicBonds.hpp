#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include "includes/MolecularMetadata/elements.hpp"

#include <string>

namespace MolecularMetadata
{
    double maxBondLengthByAtomType(Element atom1Element, Element atom2Element);
    double specificBondLength(std::string query1, std::string query2);

} // namespace MolecularMetadata

#endif
