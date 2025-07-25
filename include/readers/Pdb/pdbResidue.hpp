#ifndef INCLUDE_READERS_PDB_PDBRESIDUE_HPP
#define INCLUDE_READERS_PDB_PDBRESIDUE_HPP

#include "include/readers/Pdb/pdbAtom.hpp"
#include "include/readers/Pdb/pdbData.hpp"
#include "include/readers/Pdb/pdbResidueId.hpp"

#include <ostream>
#include <sstream>
#include <string>

namespace gmml
{
    namespace pdb
    {
        ResidueId pdbResidueId(const PdbData& data, size_t residueId);
        std::string residueStringId(const PdbData& data, size_t residueId);
        size_t readResidue(
            PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine);

        size_t readResidue(
            PdbData& data,
            size_t moleculeId,
            const std::string& name,
            ResidueType type,
            bool hasTerCard,
            size_t referenceResidue);

        std::string getNumberAndInsertionCode(const PdbData& data, size_t residueId);
    } // namespace pdb
} // namespace gmml

#endif
