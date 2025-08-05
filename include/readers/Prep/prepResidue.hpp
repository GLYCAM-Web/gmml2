#ifndef INCLUDE_READERS_PREP_PREPRESIDUE_HPP
#define INCLUDE_READERS_PREP_PREPRESIDUE_HPP

#include "include/readers/Prep/prepAtom.hpp"
#include "include/readers/Prep/prepDataTypes.hpp"

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        void initializePrepResidue(PrepData& data, std::istream& in_file, std::string& line);

        void generate3dStructure(PrepData& data, size_t index);
        void deleteDummyAtoms(PrepData& data, size_t index);
        void setConnectivities(PrepData& data, size_t index);
        std::vector<std::string> getAtomNames(PrepData& data, size_t residueId);
        std::vector<std::string> getHeavyAtomNames(PrepData& data, size_t residueId);
        double calculatePrepResidueCharge(PrepData& data, size_t residueId);

        std::string residueToString(const PrepData& data, size_t index);
        void writeResidue(const PrepData& data, size_t index, std::ostream& stream);
    } // namespace prep
} // namespace gmml

#endif
