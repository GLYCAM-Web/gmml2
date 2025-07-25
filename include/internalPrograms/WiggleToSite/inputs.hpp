#ifndef INCLUDE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
#define INCLUDE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
#include "include/readers/Pdb/pdbResidueId.hpp"

#include <string>

namespace gmml
{
    struct WiggleToSiteInputs
    {
        // ctor
        WiggleToSiteInputs(std::string inputFileName);
        // Members
        std::string carbohydrateSequence_ = "";
        uint carbohydrateSuperimpositionResidue_ = 0;
        uint carbohydrateWigglingResidue_ = 0;
        std::string substrateFile_ = "";
        pdb::ResidueId superimpositionTargetResidue_;
        pdb::ResidueId wigglingTargetResidue_;
        int substrateModelNumber_ = 1;
        int persistCycles_ = 5;
        bool isDeterministic_ = false;
        std::string Print();
    };
} // namespace gmml

#endif
