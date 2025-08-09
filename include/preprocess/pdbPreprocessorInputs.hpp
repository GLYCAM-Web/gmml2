#ifndef INCLUDE_PREPROCESS_PDBPREPROCESSORINPUTS_HPP
#define INCLUDE_PREPROCESS_PDBPREPROCESSORINPUTS_HPP

#include "include/fileType/pdb/pdbResidueId.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace preprocess
    {
        struct PreprocessorOptions
        {
            std::string chainNTermination;
            std::string chainCTermination;
            std::string gapNTermination;
            std::string gapCTermination;
            std::vector<std::pair<std::string, std::string>>
                hisSelections; // e.g. pair: residue id like this <"HIS_20_?_A_1", "HID">
        };

        static const PreprocessorOptions defaultPreprocessorOptions = {"NH3+", "CO2-", "COCH3", "NHCH3", {}};

        struct DisulfideBond
        {
            pdb::ResidueId residue1;
            pdb::ResidueId residue2;
            double distance;
        };

        struct GapInAminoAcidChain
        {
            std::string chainId;
            std::string residueBeforeGap;
            std::string residueAfterGap;
            std::string terminationBeforeGap;
            std::string terminationAfterGap;
        };

        struct AtomInfo
        {
            std::string name;
            pdb::ResidueId residue;
        };

        struct ChainTerminal
        {
            std::string chainId;
            std::string startIndex;
            std::string endIndex;
            std::string nTermination;
            std::string cTermination;
        };

        struct PreprocessorInformation
        {
            std::vector<AtomInfo> unrecognizedAtoms;
            std::vector<AtomInfo> missingHeavyAtoms;
            std::vector<pdb::ResidueId> unrecognizedResidues;
            std::vector<GapInAminoAcidChain> missingResidues;
            std::vector<ChainTerminal> chainTerminals;
            std::vector<pdb::ResidueId> hisResidues;
            std::vector<DisulfideBond> cysBondResidues;
            std::vector<pdb::ResidueId> nonNaturalProteinResidues;
        };
    } // namespace preprocess
} // namespace gmml

#endif
