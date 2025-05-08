#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_BONDBYDISTANCE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_BONDBYDISTANCE_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>

namespace cds
{
    void bondAtomsByDistance(std::vector<cds::Atom*> atoms);
    void bondAtomsAndResiduesByDistance(pdb::PdbData& pdbData, std::vector<cds::Residue*> residues);
    void distanceBondIntra(std::vector<cds::Residue*> residues);
    void distanceBondInter(pdb::PdbData& pdbData, std::vector<cds::Residue*> residues);
} // namespace cds
#endif
