#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBATOM_HPP
// See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM for an explanation of atom
// formats in PDB files

#include "includes/CentralDataStructure/atom.hpp"
#include <string>

namespace pdb
{
    void readAtom(cds::Atom* atom, const std::string& line);
} // namespace pdb
#endif
