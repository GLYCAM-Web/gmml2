#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP

#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    unsigned long long generateResidueLinkageIndex();
    ResidueLink findResidueLink(std::pair<cds::Residue*, cds::Residue*> residues);
    void determineAtomsThatMove(std::vector<RotatableDihedral>& dihedrals);
    void determineResiduesForOverlapCheck(ResidueLinkage& linkage);
    ResidueLinkage createResidueLinkage(ResidueLink& link);
} // namespace cds
#endif
