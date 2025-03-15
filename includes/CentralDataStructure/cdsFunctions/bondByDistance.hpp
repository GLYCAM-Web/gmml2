#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_BONDBYDISTANCE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_BONDBYDISTANCE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include <vector>

namespace cds
{
    void bondAtomsByDistance(std::vector<cds::Atom*> atoms);
    void bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues);
    void distanceBondIntra(std::vector<cds::Residue*> residues);
    void distanceBondInter(std::vector<cds::Residue*> residues);
} // namespace cds
#endif
