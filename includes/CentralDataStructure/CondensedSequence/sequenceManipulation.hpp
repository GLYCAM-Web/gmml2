#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    size_t instantiateNode(SequenceData& result, size_t parent, const AbstractSequence& data, size_t id);
    SequenceData instantiate(const AbstractSequence& data);
    SequenceData rearrange(const SequenceData& sequence, const std::vector<size_t>& residueOrder,
                           const std::vector<size_t>& edgeOrder);
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    SequenceData reordered(const SequenceData& sequence);
} // namespace cdsCondensedSequence
#endif
