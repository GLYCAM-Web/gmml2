#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/Assembly/assemblyTypes.hpp"

#include <ostream>
#include <string>
#include <array>
#include <vector>

namespace cds
{
    void writeAssemblyToPdb(std::ostream& stream, const assembly::Graph& graph,
                            const std::vector<std::vector<size_t>>& residueIndices,
                            const std::vector<std::vector<bool>>& residueTER,
                            const std::vector<std::array<size_t, 2>>& connectionIndices, const PdbFileData& data);
    void writeMoleculeToPdb(std::ostream& stream, const assembly::Graph& graph,
                            const std::vector<size_t>& residueIndices, const std::vector<bool>& residueTER,
                            const PdbFileData& data);
    void writeAtomToPdb(std::ostream& stream, const PdbFileFormat& format, const PdbFileResidueData& residues,
                        size_t residueIndex, const PdbFileAtomData& atoms, size_t atomIndex);
    void writeConectCards(std::ostream& stream, const std::vector<uint>& atomNumbers,
                          const std::vector<std::array<size_t, 2>>& connectionIndices);
    void theEnd(std::ostream& stream);

} // namespace cds
#endif
