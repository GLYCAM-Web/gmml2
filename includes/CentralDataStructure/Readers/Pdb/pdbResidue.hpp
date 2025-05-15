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

    size_t addResidue(PdbData& data, size_t moleculeId, size_t position, const ResidueEntry& entry);
    size_t readResidue(PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine);
    size_t readResidue(PdbData& data, size_t moleculeId, const std::string& name, cds::ResidueType type,
                       bool hasTerCard, size_t referenceResidue);
    std::string getNumberAndInsertionCode(const PdbData& data, size_t residueId);

    size_t addPdbAtom(PdbData& data, size_t residueId, const AtomEntry& entry);
    size_t addPdbAtom(PdbData& data, size_t residueId, const std::string& line);
    size_t addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& coordinate);
    void deletePdbAtom(PdbData& data, size_t residueId, size_t atomId);

    void modifyNTerminal(PdbData& data, size_t residueId, const std::string& type);
    void modifyCTerminal(PdbData& data, size_t residueId, const std::string& type);
} // namespace pdb
#endif
