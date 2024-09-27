#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <utility>
#include <vector>

namespace cds
{
    std::vector<std::pair<int, int>> atomPairNumbers(const std::vector<std::pair<Atom*, Atom*>>& pairs);
    std::vector<std::pair<size_t, size_t>> atomPairVectorIndices(const std::vector<Atom*>& atoms,
                                                                 const std::vector<std::pair<Atom*, Atom*>>& pairs);
    std::vector<std::pair<Atom*, Atom*>> atomPairsConnectedToOtherResidues(std::vector<Atom*> atoms);
    std::vector<std::pair<Atom*, Atom*>> atomPairsConnectedToOtherResidues(std::vector<Residue*> residues);
    std::vector<Atom*> atomsConnectedToOtherResidues(std::vector<Atom*> atoms);
    void setIntraConnectivity(std::vector<Residue*> residues);
    void setInterConnectivity(std::vector<Residue*> residues);
} // namespace cds

#endif
