#ifndef INCLUDE_OFF_OFFFILEWRITER_HPP
#define INCLUDE_OFF_OFFFILEWRITER_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/off/offFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace off
    {
        void WriteOffFileUnit(
            std::ostream& stream,
            const OffFileFormat& format,
            const assembly::Graph& graph,
            const OffFileResidueData& residues,
            const OffFileAtomData& atoms,
            const std::vector<size_t>& residueIndices,
            const std::string& unitName);

        void WriteResiduesIndividuallyToOffFile(
            std::ostream& stream, const assembly::Graph& graph, const OffFileData& data);

        void WriteResiduesTogetherToOffFile(
            std::ostream& stream,
            const assembly::Graph& graph,
            const OffFileData& data,
            const std::vector<size_t>& residueIndices,
            const std::string& unitName);
    } // namespace off
} // namespace gmml

#endif
