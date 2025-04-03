#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_PSIANGLEHYDROGEN_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_PSIANGLEHYDROGEN_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <vector>

namespace cds
{
    Atom* findHydrogenForPsiAngle(const Atom* atom);
    void createHydrogenForPsiAngles(const GlycamMetadata::DihedralAngleDataTable& metadataTable, Residue* residue,
                                    std::vector<DihedralAtoms>& dihedralAtoms,
                                    const std::vector<std::vector<size_t>>& metadataIndices);
} // namespace cds
#endif
