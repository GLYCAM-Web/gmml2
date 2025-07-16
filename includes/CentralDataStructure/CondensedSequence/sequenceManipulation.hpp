#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    struct ChainState
    {
        size_t tail;
        size_t head;
    };

    struct SequenceAndLinkageData
    {
        SequenceData data;
        std::vector<std::string> linkage;
    };

    ChainState instantiateNode(SequenceAndLinkageData& result, const ChainState& state, const AbstractSequence& data,
                               size_t id);
    SequenceData instantiate(const AbstractSequence& data);
    SequenceData rearrange(const SequenceData& sequence, const std::vector<size_t>& residueOrder,
                           const std::vector<size_t>& edgeOrder);
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    SequenceData reordered(const SequenceData& sequence);
} // namespace cdsCondensedSequence
#endif
