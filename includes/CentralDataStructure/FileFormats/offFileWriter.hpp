#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_OFFFILEWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    void WriteOffFileUnit(std::ostream& stream, const OffFileFormat& format, const assembly::Graph& graph,
                          const OffFileResidueData& residues, const OffFileAtomData& atoms,
                          const std::vector<size_t>& residueIndices, const std::string& unitName);
    void WriteResiduesIndividuallyToOffFile(std::ostream& stream, const assembly::Graph& graph,
                                            const OffFileData& data);
    void WriteResiduesTogetherToOffFile(std::ostream& stream, const assembly::Graph& graph, const OffFileData& data,
                                        const std::string& unitName);
} // namespace cds
#endif
