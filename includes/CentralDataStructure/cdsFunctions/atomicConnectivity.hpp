#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP

#include "includes/CentralDataStructure/residue.hpp"

#include <vector>

namespace cds
{
    void setIntraConnectivity(std::vector<cds::Residue*> residues);
    void setInterConnectivity(std::vector<cds::Residue*> residues);
} // namespace cds

#endif
