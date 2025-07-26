#ifndef INCLUDE_CENTRALDATASTRUCTURE_PRINT_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_PRINT_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/geometryTypes.hpp"

#include <ostream>

namespace gmml
{
    void print(std::ostream& out, const Coordinate& coord);
    void print(std::ostream& out, const Atom& atom);
    void print(std::ostream& out, const Residue& residue);
} // namespace gmml
#endif
