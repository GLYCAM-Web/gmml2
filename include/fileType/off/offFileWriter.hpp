#ifndef INCLUDE_FILETYPE_OFF_OFFFILEWRITER_HPP
#define INCLUDE_FILETYPE_OFF_OFFFILEWRITER_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/off/offFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace off
    {
        void writeResiduesIndividually(std::ostream& stream, const assembly::Graph& graph, const OffFileData& data);

        void writeResiduesTogether(
            std::ostream& stream,
            const assembly::Graph& graph,
            const OffFileData& data,
            const std::vector<size_t>& residueIndices,
            const std::string& unitName);
    } // namespace off
} // namespace gmml

#endif
