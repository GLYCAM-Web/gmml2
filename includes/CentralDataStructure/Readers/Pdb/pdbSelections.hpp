#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBSELECTIONS_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"

namespace pdb
{
    std::vector<cds::Atom*> getAtoms(const std::vector<PdbModel>& models);
    std::vector<cds::Residue*> getResidues(const std::vector<PdbModel>& models);
    PdbResidue* residueSelector(const PdbFile& pdbFile, const pdb::ResidueId& residueId, const int modelNumber = 0);
    PdbResidue* residueSelector(std::vector<cds::Residue*> residues, const pdb::ResidueId& queryId);
} // namespace pdb
#endif
