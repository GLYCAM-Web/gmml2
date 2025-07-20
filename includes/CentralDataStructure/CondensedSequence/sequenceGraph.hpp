#ifndef INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP
#define INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    std::vector<std::string> sequenceDerivatives(const SequenceData& sequence);
    std::vector<std::string> sequenceMonosaccharideNames(const SequenceData& sequence);
    graph::Graph condensedSequenceGraph(const SequenceData& sequence);
} // namespace cdsCondensedSequence

#endif
