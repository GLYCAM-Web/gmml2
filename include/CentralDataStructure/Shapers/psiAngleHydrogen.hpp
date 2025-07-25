#ifndef INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_PSIANGLEHYDROGEN_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_PSIANGLEHYDROGEN_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <vector>

namespace gmml
{
    Atom* findHydrogenForPsiAngle(const Atom* atom);
    void createHydrogenForPsiAngles(
        const DihedralAngleDataTable& metadataTable,
        Residue* residue,
        std::vector<DihedralAtoms>& dihedralAtoms,
        const std::vector<std::vector<size_t>>& metadataIndices);
} // namespace gmml
#endif
