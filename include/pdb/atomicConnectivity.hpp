#ifndef INCLUDE_PDB_ATOMICCONNECTIVITY_HPP
#define INCLUDE_PDB_ATOMICCONNECTIVITY_HPP

#include "include/metadata/aminoAcids.hpp"
#include "include/pdb/pdbData.hpp"

namespace gmml
{
    namespace pdb
    {
        void setIntraConnectivity(const AminoAcidTable& table, PdbData& data);
        void setInterConnectivity(const AminoAcidTable& table, PdbData& pdbData);
    } // namespace pdb
} // namespace gmml

#endif
