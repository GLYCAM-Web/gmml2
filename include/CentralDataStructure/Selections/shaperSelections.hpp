#ifndef INCLUDE_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_SELECTIONS_SHAPERSELECTIONS_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    bool FindPathBetweenTwoAtoms(
        Atom* current_atom,
        Residue* currentResidue,
        Atom* target_atom,
        Residue* targetResidue,
        std::vector<Atom*>* atom_path,
        bool* found);
    void FindAtomsConnectingResidues(
        Atom* current_atom,
        const Residue* currentResidue,
        const Residue* otherResidue,
        std::vector<Atom*>* connecting_atoms,
        bool* found_neighbor);
    void ClearAtomLabels(Residue* residue);
    ResidueLinkage* selectLinkageWithIndex(
        std::vector<ResidueLinkage>& inputLinkages, const long long unsigned int indexQuery);
    std::vector<DihedralAtoms> splitAtomVectorIntoRotatableDihedrals(const std::vector<Atom*>& atoms);
    std::vector<DihedralAtoms> findRotatableDihedralsConnectingResidues(const ResidueLink& link);
} // namespace gmml

#endif
