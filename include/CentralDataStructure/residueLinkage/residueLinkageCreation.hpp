#ifndef INCLUDE_CENTRALDATASTRUCTURE_RESIDUELINKAGE_RESIDUELINKAGECREATION_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_RESIDUELINKAGE_RESIDUELINKAGECREATION_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <utility>
#include <vector>

namespace gmml
{
    unsigned long long generateResidueLinkageIndex();
    ResidueLink findResidueLink(std::pair<Residue*, Residue*> residues);
    void determineAtomsThatMove(std::vector<RotatableBond>& bonds);
    ResidueLinkage createResidueLinkage(const DihedralAngleDataTable& metadataTable, ResidueLink& link);
} // namespace gmml
#endif
