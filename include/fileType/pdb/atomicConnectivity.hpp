#ifndef INCLUDE_FILETYPE_PDB_ATOMICCONNECTIVITY_HPP
#define INCLUDE_FILETYPE_PDB_ATOMICCONNECTIVITY_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/metadata/aminoAcids.hpp"

namespace gmml
{
    namespace pdb
    {
        void setIntraConnectivity(const AminoAcidTable& table, PdbData& data);
        void setInterConnectivity(const AminoAcidTable& table, PdbData& pdbData, const assembly::Bounds& bounds);
    } // namespace pdb
} // namespace gmml

#endif
