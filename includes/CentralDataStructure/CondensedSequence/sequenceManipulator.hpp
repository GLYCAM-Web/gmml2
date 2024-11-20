#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"

#include <memory>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    std::vector<ParsedResidue*> parsedResiduesOrderedByConnectivity(std::vector<ParsedResidue*> residues);
    void setIndexByConnectivity(std::vector<cdsCondensedSequence::ParsedResidue*> residues);
    void labelSequence(std::vector<ParsedResidue*> residues);
    std::string printGraphViz(GraphVizDotConfig& configs, std::vector<ParsedResidue*> residues);
    std::string printSequence(std::vector<ParsedResidue*> residues, bool withLabels = false,
                              bool iupacCondensed = false);
    std::string reorderSequence(std::vector<std::unique_ptr<ParsedResidue>>& residues);
    ParsedResidue* terminalResidue(std::vector<ParsedResidue*> residues);
} // namespace cdsCondensedSequence
#endif
