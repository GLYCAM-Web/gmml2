#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_ATOMICCONNECTIVITY_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"

#include <utility>
#include <array>
#include <vector>

namespace cds
{
    std::vector<std::pair<int, int>> atomPairNumbers(const std::vector<std::pair<Atom*, Atom*>>& pairs);
    std::vector<std::array<size_t, 2>> atomPairVectorIndices(const std::vector<Atom*>& atoms,
                                                             const std::vector<std::pair<Atom*, Atom*>>& pairs);
    std::vector<std::pair<Atom*, Atom*>> atomPairsConnectedToOtherResidues(std::vector<Atom*> atoms);
    std::vector<std::pair<Atom*, Atom*>> atomPairsConnectedToOtherResidues(std::vector<Residue*> residues);
    std::vector<Atom*> atomsConnectedToOtherResidues(std::vector<Atom*> atoms);
    void setIntraConnectivity(const MolecularMetadata::AminoAcidTable& table, std::vector<Residue*> residues);
    void setInterConnectivity(const MolecularMetadata::AminoAcidTable& table, pdb::PdbData& pdbData,
                              std::vector<Residue*> residues);
} // namespace cds

#endif
