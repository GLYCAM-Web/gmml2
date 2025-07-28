#ifndef INCLUDE_PDB_PDBATOM_HPP
#define INCLUDE_PDB_PDBATOM_HPP
// See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM for an explanation of atom
// formats in PDB files

#include "include/pdb/pdbData.hpp"

#include <string>

namespace gmml
{
    namespace pdb
    {
        AtomEntry readAtom(const std::string& line);
    } // namespace pdb
} // namespace gmml

#endif
