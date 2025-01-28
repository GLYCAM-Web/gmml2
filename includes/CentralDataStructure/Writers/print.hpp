#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PRINT_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PRINT_HPP

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <ostream>

namespace cds
{
    void print(std::ostream& out, const cds::Coordinate& coord);
    void print(std::ostream& out, const cds::Atom& atom);
    void print(std::ostream& out, const cds::Residue& residue);
} // namespace cds
#endif
