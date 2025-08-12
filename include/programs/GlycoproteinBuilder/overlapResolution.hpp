#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_OVERLAPRESOLUTION_HPP

#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/metadata/sidechainRotamers.hpp"
#include "include/programs/GlycoproteinBuilder/gpBuilderStats.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct ResolutionSettings
        {
            uint64_t rngSeed;
            ulong persistCycles;
            ulong numberOfSamples;
            bool useInitialGlycositeResidueConformation;
            bool moveOverlappingSidechains;
            bool deleteSitesUntilResolved;
            bool prepareForMD;
            bool writeOffFile;
        };

        void resolveOverlaps(
            std::vector<StructureStats>& stats,
            const SidechainRotamerData& sidechainRotamers,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const GlycoproteinAssembly& assembly,
            const OverlapSettings& overlapSettings,
            const ResolutionSettings& resolutionSettings,
            const std::string& outputDir,
            const std::vector<std::string>& headerLines,
            int numThreads);
    } // namespace gpbuilder
} // namespace gmml

#endif
