#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEUTIL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEUTIL_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/residueTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    size_t parentEdge(const SequenceData& sequence, size_t nodeId);
    std::string edgeLinkage(const SequenceData& sequence, size_t edgeId);
    std::string mainLinkage(cds::ResidueType type, const std::string& linkage);
    std::string mainLinkage(const SequenceData& sequence, size_t residueId);
} // namespace cdsCondensedSequence
#endif
