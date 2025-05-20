#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBFUNCTIONS_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"

#include <string>

namespace pdb
{
    int checkShiftFromSerialNumberOverrun(const std::string& line);
    int checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift = 0);
    cds::Coordinate checkShiftsAndExtractCoordinate(const std::string& line);
    cds::Coordinate coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z);
    void expandLine(std::string& line, int length);
    size_t addAtom(PdbData& data, size_t residueId, const AtomEntry& entry);
    size_t addAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& coordinate);
    void deleteAtom(PdbData& data, size_t residueId, size_t atomId);
    size_t addResidue(PdbData& data, size_t moleculeId, size_t position, const ResidueEntry& entry);
    void addBond(PdbData& data, size_t atom1, size_t atom2);
    size_t findResidueAtom(const PdbData& data, size_t residueId, const std::string& atomName);
}; // namespace pdb
#endif
