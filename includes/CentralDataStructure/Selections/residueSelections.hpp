#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_RESIDUESELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_RESIDUESELECTIONS_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"

#include <vector>

namespace cdsSelections
{
    using cds::Residue;
    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, cds::ResidueType queryType, const bool invert = false);
    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, std::vector<cds::ResidueType> queryTypes, const bool invert = false);
    void FindConnectedResidues(std::vector<Residue*>& visitedList, Residue* current);
} // namespace cdsSelections
#endif
