#ifndef INCLUDE_SEQUENCE_SEQUENCEUTIL_HPP
#define INCLUDE_SEQUENCE_SEQUENCEUTIL_HPP

#include "include/CentralDataStructure/residueTypes.hpp"
#include "include/sequence/sequenceTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        size_t parentEdge(const SequenceData& sequence, size_t nodeId);
        std::string edgeLinkage(const SequenceData& sequence, size_t edgeId);
        std::string mainLinkage(ResidueType type, const std::string& linkage);
        std::string mainLinkage(const SequenceData& sequence, size_t residueId);
    } // namespace sequence
} // namespace gmml

#endif
