#ifndef INCLUDE_PDB_PDBSELECTIONS_HPP
#define INCLUDE_PDB_PDBSELECTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbResidueId.hpp"

namespace gmml
{
    namespace pdb
    {
        std::vector<Atom*> getAtoms(const std::vector<Assembly*>& assemblies);
        std::vector<Residue*> getResidues(const std::vector<Assembly*>& assemblies);
        size_t residueSelector(const PdbData& data, const pdb::ResidueId& residueId, const int modelNumber = 0);
        size_t residueSelector(const PdbData& data, std::vector<size_t> residueIds, const pdb::ResidueId& queryId);
    } // namespace pdb
} // namespace gmml

#endif
