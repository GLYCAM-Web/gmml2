#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPRINTER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPRINTER_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <functional>
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

    SequencePrintConfig defaultConfig();
    SequencePrintConfig defaultConfigUnsorted();
    SequencePrintConfig defaultConfigLabelled();
    SequencePrintConfig iupacConfig();
    std::string printSequence(const SequencePrintConfig& config, const SequenceData& sequence);
    std::string printGraphViz(GraphVizDotConfig& configs, const SequenceData& sequence);
} // namespace cdsCondensedSequence
#endif
