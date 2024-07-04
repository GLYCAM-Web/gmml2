#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"

namespace pdb
{
    PdbResidue* residueSelector(const PdbFile& pdbFile, const pdb::ResidueId& residueId, const int modelNumber = 0);
    PdbResidue* residueSelector(std::vector<cds::Residue*> residues, const pdb::ResidueId& queryId);
} // namespace pdb
#endif
