#ifndef INCLUDE_PDB_PDBFUNCTIONS_HPP
#define INCLUDE_PDB_PDBFUNCTIONS_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/pdb/pdbData.hpp"

#include <string>

namespace gmml
{
    namespace pdb
    {
        std::string residueParameterName(const PdbData& data, size_t residueId);
        int checkShiftFromSerialNumberOverrun(const std::string& line);
        int checkSecondShiftFromResidueNumberOverrun(const std::string& line, const int shift = 0);
        Coordinate checkShiftsAndExtractCoordinate(const std::string& line);
        Coordinate coordinateFromStrings(const std::string& x, const std::string& y, const std::string& z);
        void expandLine(std::string& line, int length);
        size_t addAtom(PdbData& data, size_t residueId, const AtomEntry& entry);
        size_t addAtom(PdbData& data, size_t residueId, const std::string& name, const Coordinate& coordinate);
        void deleteAtom(PdbData& data, size_t residueId, size_t atomId);
        size_t addResidue(PdbData& data, size_t moleculeId, size_t position, const ResidueEntry& entry);
        void addBond(PdbData& data, size_t atom1, size_t atom2);
        size_t findResidueAtom(const PdbData& data, size_t residueId, const std::string& atomName);
    }; // namespace pdb
} // namespace gmml

#endif
