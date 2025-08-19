#ifndef INCLUDE_METADATA_ATOMICBONDS_HPP
#define INCLUDE_METADATA_ATOMICBONDS_HPP

#include "include/metadata/elements.hpp"

#include <string>

namespace gmml
{
    double deoxyHydrogenBondLength();
    double hydrogenCovalentBondMaxLength();
    double maxBondLengthByAtomType(Element atom1Element, Element atom2Element);
    double specificBondLength(const std::string& query1, const std::string& query2);
} // namespace gmml

#endif
