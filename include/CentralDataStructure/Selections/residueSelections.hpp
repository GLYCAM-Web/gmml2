#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_RESIDUESELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_RESIDUESELECTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/readers/Pdb/pdbData.hpp"

#include <vector>

namespace gmml
{
    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, ResidueType queryType, const bool invert = false);
    std::vector<Residue*> selectResiduesByType(
        std::vector<Residue*> inputResidues, std::vector<ResidueType> queryTypes, const bool invert = false);
    void FindConnectedResidues(std::vector<Residue*>& visitedList, Residue* current);
} // namespace gmml

#endif
