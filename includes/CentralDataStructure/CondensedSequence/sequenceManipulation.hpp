#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATION_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    std::string mainLinkage(cds::ResidueType type, const std::string& linkage);
    std::string mainLinkage(const SequenceData& sequence, size_t residueId);
    SequenceData instantiate(const AbstractSequence& data);
    std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds);
    SequenceData reordered(const SequenceData& sequence);
} // namespace cdsCondensedSequence
#endif
