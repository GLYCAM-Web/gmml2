#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"

#include <string>
#include <sstream>
#include <ostream>

namespace pdb
{
    ResidueId pdbResidueId(const PdbData& data, size_t residueId);
    std::string residueStringId(const PdbData& data, size_t residueId);

    size_t readResidue(PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine);
    size_t readResidue(PdbData& data, size_t moleculeId, const std::string& name, cds::ResidueType type,
                       bool hasTerCard, size_t referenceResidue);
    std::string getNumberAndInsertionCode(const PdbData& data, size_t residueId);
} // namespace pdb
#endif
