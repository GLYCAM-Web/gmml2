#ifndef INCLUDES_READERS_PDB_ATOMICCONNECTIVITY_HPP
#define INCLUDES_READERS_PDB_ATOMICCONNECTIVITY_HPP

#include "include/metadata/aminoAcids.hpp"
#include "include/readers/Pdb/pdbData.hpp"

#include <array>
#include <utility>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        void setIntraConnectivity(const AminoAcidTable& table, PdbData& data);
        void setInterConnectivity(const AminoAcidTable& table, PdbData& pdbData);
    } // namespace pdb
} // namespace gmml

#endif
