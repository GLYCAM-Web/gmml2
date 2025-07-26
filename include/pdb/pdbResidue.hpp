#ifndef INCLUDE_PDB_PDBRESIDUE_HPP
#define INCLUDE_PDB_PDBRESIDUE_HPP

#include "include/pdb/pdbAtom.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbResidueId.hpp"

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
