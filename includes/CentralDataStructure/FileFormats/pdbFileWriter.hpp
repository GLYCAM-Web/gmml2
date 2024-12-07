#ifndef INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_FILEFORMATS_PDBFILEWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace cds
{
    void writeMoleculeToPdb(std::ostream& stream, const std::vector<size_t>& residueIndices,
                            const std::vector<bool>& residueTER, const PdbFileData& data);
    void writeAtomToPdb(std::ostream& stream, const PdbFileResidueData& residues, size_t residueIndex,
                        const PdbFileAtomData& atoms, size_t atomIndex);
    void writeConectCards(std::ostream& stream, const std::vector<int>& atomNumbers,
                          std::vector<std::pair<size_t, size_t>> connectionIndices);

} // namespace cds
#endif
