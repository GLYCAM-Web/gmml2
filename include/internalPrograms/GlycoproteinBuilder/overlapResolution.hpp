#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP

#include "include/geometry/overlap.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/metadata/sidechainRotamers.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        void resolveOverlaps(
            const SidechainRotamerData& sidechainRotamers,
            const GlycoproteinAssembly& assembly,
            const GlycoproteinBuilderInputs& settings,
            uint64_t rngSeed,
            const std::string& outputDir,
            const std::vector<std::string>& headerLines,
            int numThreads);
    }
} // namespace gmml

#endif
