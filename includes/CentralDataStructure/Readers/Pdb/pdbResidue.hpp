#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
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

    void readResidue(PdbData& data, size_t residueId, std::stringstream& singleResidueSecion, std::string firstLine);
    void readResidue(PdbData& data, size_t residueId, size_t referenceResidue);
    std::string getNumberAndInsertionCode(const PdbData& data, size_t residueId);

    size_t addPdbAtom(PdbData& data, size_t residueId, const std::string& line);
    size_t addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& c);
    void deletePdbAtom(PdbData& data, size_t residueId, size_t atomId);

    void modifyNTerminal(PdbData& data, size_t residueId, const std::string& type);
    void modifyCTerminal(PdbData& data, size_t residueId, const std::string& type);
} // namespace pdb
#endif
