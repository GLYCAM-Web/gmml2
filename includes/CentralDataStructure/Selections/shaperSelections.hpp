#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP

#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include <vector>
#include <string>

namespace cdsSelections
{
    bool FindPathBetweenTwoAtoms(cds::Atom* current_atom, cds::Residue* currentResidue, cds::Atom* target_atom,
                                 cds::Residue* targetResidue, std::vector<cds::Atom*>* atom_path, bool* found);
    void FindAtomsConnectingResidues(cds::Atom* current_atom, const cds::Residue* currentResidue,
                                     const cds::Residue* otherResidue, std::vector<cds::Atom*>* connecting_atoms,
                                     bool* found_neighbor);
    void ClearAtomLabels(cds::Residue* residue);
    cds::ResidueLinkage* selectLinkageWithIndex(std::vector<cds::ResidueLinkage>& inputLinkages,
                                                const long long unsigned int indexQuery);
    std::vector<cds::ResidueLinkage> SplitLinkagesIntoPermutants(std::vector<cds::ResidueLinkage>& inputLinkages);
    std::vector<cds::DihedralAtoms> splitAtomVectorIntoRotatableDihedrals(const std::vector<cds::Atom*>& atoms);
    std::vector<cds::DihedralAtoms> findRotatableDihedralsConnectingResidues(const cds::ResidueLink& link);
} // namespace cdsSelections
#endif
