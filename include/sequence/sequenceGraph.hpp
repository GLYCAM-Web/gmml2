#ifndef INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP
#define INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/sequence/sequenceTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        std::vector<std::string> sequenceDerivatives(const SequenceData& sequence);
        std::vector<std::string> sequenceMonosaccharideNames(const SequenceData& sequence);
        graph::Graph condensedSequenceGraph(const SequenceData& sequence);
    } // namespace sequence
} // namespace gmml

#endif
