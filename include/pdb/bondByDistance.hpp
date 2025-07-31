#ifndef INCLUDE_PDB_BONDBYDISTANCE_HPP
#define INCLUDE_PDB_BONDBYDISTANCE_HPP

#include "include/pdb/pdbData.hpp"

#include <vector>

namespace gmml
{
    namespace pdb
    {
        bool isWithinBondingDistance(const PdbData& data, size_t atomA, size_t atomB);
        void bondAtomsByDistance(PdbData& data, const std::vector<size_t>& atoms);
        void bondAtomsAndResiduesByDistance(PdbData& data, const assembly::Bounds& bounds);
        void distanceBondIntra(PdbData& data, const std::vector<size_t>& residues);
        void distanceBondInter(PdbData& data, const assembly::Bounds& bounds, const std::vector<size_t>& residues);
    } // namespace pdb
} // namespace gmml

#endif
