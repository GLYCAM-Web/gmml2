#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <vector>

namespace cdsCondensedSequence
{
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    SequenceData pruned(const SequenceData& sequence);
    SequenceData reordered(const SequenceData& sequence);
    SequenceData constructSequence(pcg32& rng, SequenceData data);

    inline SequenceData parseAndReorder(const std::string& sequence)
    {
        pcg32 rng(0);
        return reordered(constructSequence(rng, parseSequence(sequence)));
    }
} // namespace cdsCondensedSequence
#endif
