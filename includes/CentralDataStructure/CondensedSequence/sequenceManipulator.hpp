#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"

#include <memory>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    std::vector<ParsedResidue*> parsedResiduesOrderedByConnectivity(std::vector<ParsedResidue*> residues);
    void setIndexByConnectivity(std::vector<cdsCondensedSequence::ParsedResidue*> residues);
    SequenceData reordered(const SequenceData& sequence);
    void createParsedResidues(std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs,
                              const SequenceData& sequence);
    void sortResidueEdges(std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs);
    ParsedResidue* terminalResidue(std::vector<ParsedResidue*> residues);
} // namespace cdsCondensedSequence
#endif
