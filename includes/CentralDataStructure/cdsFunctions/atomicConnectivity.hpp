#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP

#include "includes/CentralDataStructure/residue.hpp"

namespace cds
{
    void setBondingForAminoAcid(cds::Residue* proteinRes);
    bool autoConnectSuccessiveResidues(cds::Residue* cTermRes, cds::Residue* nTermRes);
    void setProteinIntraConnectivity(std::vector<cds::Residue*> proteinResidues);
    void setProteinInterConnectivity(std::vector<cds::Residue*> proteinResidues);
    void setIntraConnectivity(std::vector<cds::Residue*> residues);
    void setInterConnectivity(std::vector<cds::Residue*> residues);
} // namespace cds

#endif
