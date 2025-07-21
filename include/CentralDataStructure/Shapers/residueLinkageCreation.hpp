#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <utility>
#include <vector>

namespace gmml
{
    unsigned long long generateResidueLinkageIndex();
    ResidueLink findResidueLink(std::pair<Residue*, Residue*> residues);
    void determineAtomsThatMove(std::vector<RotatableDihedral>& dihedrals);
    ResidueLinkage createResidueLinkage(const DihedralAngleDataTable& metadataTable, ResidueLink& link);
} // namespace gmml
#endif
