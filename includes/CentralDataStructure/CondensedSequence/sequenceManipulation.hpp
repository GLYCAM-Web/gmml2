#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    SequenceData reordered(const SequenceData& sequence);

    inline SequenceData parseAndReorder(const std::string& sequence)
    {
        return reordered(parseSequence(sequence));
    }
} // namespace cdsCondensedSequence
#endif
