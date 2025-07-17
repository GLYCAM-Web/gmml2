#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGECREATION_HPP

#include "includes/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <vector>
#include <utility>

namespace cds
{
    unsigned long long generateResidueLinkageIndex();
    ResidueLink findResidueLink(std::pair<Residue*, Residue*> residues);
    void determineAtomsThatMove(std::vector<RotatableDihedral>& dihedrals);
    ResidueLinkage createResidueLinkage(const GlycamMetadata::DihedralAngleDataTable& metadataTable, ResidueLink& link);
} // namespace cds
#endif
