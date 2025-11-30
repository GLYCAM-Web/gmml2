#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP

#include "include/util/constants.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        static const std::string proteinParameter = "Protein";
        static const std::string numberOfSamplesParameter = "numberOfSamples";
        static const std::string persistCyclesParameter = "persistCycles";
        static const std::string useInitialGlycositeResidueConformationParameter =
            "useInitialGlycositeResidueConformation";
        static const std::string moveOverlappingSidechainsParameter = "moveOverlappingSidechains";
        static const std::string deleteIncompatibleSitesParameter = "deleteUnresolvableGlycosites";
        static const std::string overlapRejectionThresholdParameter = "overlapRejectionThreshold";
        static const std::string ignoreHydrogenParameter = "ignoreHydrogen";
        static const std::string seedParameter = "rngSeed";
        static const std::string prepareForMDParameter = "prepareForMD";

        struct GlycositeInput
        {
            std::string proteinResidueId = ""; // E.g. ?_20_A if no chain ID and residue number is 20 and insertion code
                                               // is A. C_20_? if chain id is C and there is no insertion code.
            std::string glycanInput = "";      // E.g. DGlcpNAcb1-4DGlcpNAcb1-OH
        };

        struct GlycoproteinBuilderInputs
        {
            std::string inputFileName;
            std::string substrateFileName = "Undefined"; // Program should throw if left as "Undefined".
            ulong numberOfSamples = 1;
            ulong persistCycles = 5;
            bool ignoreHydrogen = false;
            bool useInitialGlycositeResidueConformation = false;
            bool moveOverlappingSidechains = false;
            bool deleteSitesUntilResolved = false;
            double overlapRejectionThreshold = 0.0;
            bool isDeterministic = false;
            uint64_t seed = 0;
            bool MDprep = false;
            bool writeOffFile = false;
            std::vector<GlycositeInput> glycositesInputVector; // No default, program will throw if uninitialized.
        };

        GlycoproteinBuilderInputs readGPInputFile(std::string inputFileName);
    } // namespace gpbuilder
} // namespace gmml

#endif
