#ifndef INCLUDE_READERS_PREP_PREPATOM_HPP
#define INCLUDE_READERS_PREP_PREPATOM_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/readers/Prep/prepDataTypes.hpp"

#include <istream>
#include <ostream>
#include <string>

namespace gmml
{
    namespace prep
    {
        void initializePrepAtom(PrepData& data, size_t residueId, const std::string& line);
        Coordinate determine3dCoordinate(PrepData& data, size_t index);
        void printAtom(const PrepData& data, size_t index, std::ostream& out);
        void writeAtom(const PrepData& data, size_t index, std::ostream& stream);
    } // namespace prep
} // namespace gmml

#endif
