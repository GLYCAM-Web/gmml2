#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"

namespace pdb
{
    std::vector<cds::Atom*> getAtoms(const std::vector<cds::Assembly*>& assemblies);
    std::vector<cds::Residue*> getResidues(const std::vector<cds::Assembly*>& assemblies);
    size_t residueSelector(const PdbData& data, const pdb::ResidueId& residueId, const int modelNumber = 0);
    size_t residueSelector(const PdbData& data, std::vector<size_t> residueIds, const pdb::ResidueId& queryId);
} // namespace pdb
#endif
