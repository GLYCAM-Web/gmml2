#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP

#include "includes/CodeUtils/constants.hpp"

#include <string>
#include <vector>

namespace glycoproteinBuilder
{
    static const std::string proteinParameter                            = "Protein";
    static const std::string numberOfStructuresParameter                 = "NumberOfOutputStructures";
    static const std::string persistCyclesParameter                      = "persistCycles";
    static const std::string freezeGlycositeResidueConformationParameter = "freezeGlycositeResidueConformation";
    static const std::string allowSidechainAdjustmentParameter           = "allowSidechainAdjustment";
    static const std::string deleteIncompatibleSitesParameter            = "deleteIncompatibleSites";
    static const std::string atomPotentialRejectionThresholdParameter    = "atomPotentialForceRejectionThreshold";
    static const std::string seedParameter                               = "seed";
    static const std::string prepareForMDParameter                       = "prepareForMD";

    struct GlycositeInput
    {
        std::string proteinResidueId = ""; // E.g. ?_20_A if no chain ID and residue number is 20 and insertion code is
                                           // A. C_20_? if chain id is C and there is no insertion code.
        std::string glycanInput      = ""; // E.g. DGlcpNAcb1-4DGlcpNAcb1-OH
    };

    struct GlycoproteinBuilderInputs
    {
        std::string inputFileName;
        std::string substrateFileName           = "Undefined"; // Program should throw if left as "Undefined".
        ulong number3DStructures                = 1;
        ulong persistCycles                     = 5;
        double overlapTolerance                 = constants::overlapTolerance;
        bool freezeGlycositeResidueConformation = false;
        bool allowSidechainAdjustment           = false;
        bool deleteSitesUntilResolved           = false;
        bool rejectExcessiveGlycanOverlaps      = false;
        double atomPotentialRejectionThreshold  = 0.0;
        bool isDeterministic                    = false;
        uint64_t seed                           = 0;
        bool MDprep                             = false;
        bool writeOffFile                       = false;
        std::vector<GlycositeInput> glycositesInputVector; // No default, program will throw if uninitialized.
    };

    GlycoproteinBuilderInputs readGPInputFile(std::string inputFileName);
} // namespace glycoproteinBuilder
#endif
