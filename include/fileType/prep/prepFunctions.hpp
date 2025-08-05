#ifndef INCLUDE_FILETYPE_PREP_PREPFUNCTIONS_HPP
#define INCLUDE_FILETYPE_PREP_PREPFUNCTIONS_HPP

#include "include/fileType/prep/prepDataTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        std::vector<size_t> residueAtoms(const PrepData& data, size_t residueId);
        std::vector<size_t> residueEdges(const PrepData& data, size_t residueId);
        std::vector<size_t> atomNeighbors(const PrepData& data, size_t atomId);
        size_t findAtom(const PrepData& data, size_t residueId, const std::string& name);
        size_t addAtom(prep::PrepData& data, size_t residueId);
        size_t copyResidue(prep::PrepData& data, const prep::PrepData& reference, size_t residue);

        inline size_t copyResidue(prep::PrepData& data, size_t residue) { return copyResidue(data, data, residue); }
    } // namespace prep
} // namespace gmml

#endif
