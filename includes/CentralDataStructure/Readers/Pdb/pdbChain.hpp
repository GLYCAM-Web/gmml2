#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <string>
#include <sstream>

// #include <functional>
namespace pdb
{
    void readChain(PdbData& data, size_t moleculeId, cds::Molecule* molecule, std::stringstream& stream_block);
    std::stringstream extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line);
    void tagTerminalResidues(cds::Molecule* molecule);
    PdbResidue* getNTerminal(cds::Molecule* molecule);
    PdbResidue* getCTerminal(cds::Molecule* molecule);
    void ModifyTerminal(PdbData& data, size_t residueId, const std::string& type);
    void InsertCap(PdbData& data, size_t moleculeId, cds::Molecule* molecule, PdbResidue& refResidue,
                   const std::string& type);
} // namespace pdb
#endif
