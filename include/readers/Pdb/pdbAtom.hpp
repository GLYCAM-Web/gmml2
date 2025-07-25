#ifndef INCLUDE_READERS_PDB_PDBATOM_HPP
#define INCLUDE_READERS_PDB_PDBATOM_HPP
// See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM for an explanation of atom
// formats in PDB files

#include "include/geometry/geometryTypes.hpp"
#include "include/readers/Pdb/pdbData.hpp"

#include <string>

namespace gmml
{
    namespace pdb
    {
        AtomEntry readAtom(const std::string& line);
    } // namespace pdb
} // namespace gmml

#endif
