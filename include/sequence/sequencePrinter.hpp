#ifndef INCLUDES_CONDENSEDSEQUENCE_SEQUENCEPRINTER_HPP
#define INCLUDES_CONDENSEDSEQUENCE_SEQUENCEPRINTER_HPP

#include "include/sequence/graphViz.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/util/strings.hpp"

#include <functional>
#include <string>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        typedef std::function<std::vector<size_t>(const SequenceData&, const std::vector<size_t>&)> EdgeOrder;
        typedef std::function<std::vector<std::string>(const SequenceData&)> ResidueDisplayString;
        typedef std::function<std::vector<util::Brackets>(const SequenceData&)> ResidueLinkageBrackets;
        typedef std::function<std::vector<std::string>(const SequenceData&, const std::vector<std::vector<size_t>>&)>
            EdgeDisplayString;

        struct SequencePrintConfig
        {
            EdgeOrder edgeOrder;
            util::Brackets derivativeBrackets;
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
    } // namespace sequence
} // namespace gmml

#endif
