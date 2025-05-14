#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_BONDBYDISTANCE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_BONDBYDISTANCE_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"

#include <vector>

namespace pdb
{
    bool isWithinBondingDistance(const PdbData& data, size_t atomA, size_t atomB);
    void bondAtomsByDistance(PdbData& data, const std::vector<size_t>& atoms);
    void bondAtomsAndResiduesByDistance(PdbData& data);
    void distanceBondIntra(PdbData& data, const std::vector<size_t>& residues);
    void distanceBondInter(PdbData& data, const std::vector<size_t>& residues);
} // namespace pdb
#endif
