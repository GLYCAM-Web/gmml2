#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <functional>
#include <memory>
#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    typedef std::function<std::vector<size_t>(const SequenceData&, const std::vector<size_t>&)> EdgeOrder;
    typedef std::function<std::vector<std::string>(const SequenceData&)> ResidueDisplayString;
    typedef std::function<std::vector<codeUtils::Brackets>(const SequenceData&)> ResidueLinkageBrackets;
    typedef std::function<std::vector<std::string>(const SequenceData&, const std::vector<std::vector<size_t>>&)>
        EdgeDisplayString;

    struct SequencePrintConfig
    {
        EdgeOrder edgeOrder;
        codeUtils::Brackets derivativeBrackets;
        std::string derivativeSeparator;
        ResidueDisplayString residueNames;
        ResidueDisplayString residueLabels;
        ResidueDisplayString parentlessResidueLabels;
        ResidueLinkageBrackets residueLinkageBrackets;
        EdgeDisplayString edgeNames;
        EdgeDisplayString edgeLabels;
    };

    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    std::vector<ParsedResidue*> parsedResiduesOrderedByConnectivity(std::vector<ParsedResidue*> residues);
    void setIndexByConnectivity(std::vector<cdsCondensedSequence::ParsedResidue*> residues);
    std::string printGraphViz(GraphVizDotConfig& configs, const SequenceData& sequence);
    SequencePrintConfig defaultConfig();
    SequencePrintConfig defaultConfigUnsorted();
    SequencePrintConfig defaultConfigLabelled();
    SequencePrintConfig iupacConfig();
    std::string printSequence(const SequenceData& sequence, const SequencePrintConfig& config);
    SequenceData reordered(const SequenceData& sequence);
    void createParsedResidues(std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs,
                              const SequenceData& sequence);
    void sortResidueEdges(std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs);
    ParsedResidue* terminalResidue(std::vector<ParsedResidue*> residues);
} // namespace cdsCondensedSequence
#endif
