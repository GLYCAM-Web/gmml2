#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_ATOMICCONNECTIVITY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_ATOMICCONNECTIVITY_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"

#include <array>
#include <utility>
#include <vector>

namespace pdb
{
    void setIntraConnectivity(const MolecularMetadata::AminoAcidTable& table, PdbData& data);
    void setInterConnectivity(const MolecularMetadata::AminoAcidTable& table, PdbData& pdbData);
} // namespace pdb

#endif
