#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"

#include <string>
#include <vector>

namespace glycoproteinBuilder
{
    void resolveOverlaps(const MolecularMetadata::SidechainRotamerData& sidechainRotamers,
                         const GlycoproteinAssembly& assembly, const GlycoproteinBuilderInputs& settings,
                         const std::string& outputDir, const std::vector<std::string>& headerLines, int numThreads);
} // namespace glycoproteinBuilder

#endif
