#ifndef INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP
#define INCLUDES_SEQUENCE_SEQUENCEGRAPH_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"

#include <vector>

namespace cdsCondensedSequence
{
    struct SequenceGraph
    {
        std::vector<cds::ResidueType> types;
        std::vector<std::string> names;
        std::vector<std::string> monosaccharideNames;
        std::vector<std::string> ringTypes;
        std::vector<std::string> linkageNames;
        std::vector<std::string> configurations;
        std::vector<std::string> linkages;
        std::vector<std::string> derivatives;
        graph::Graph graph;
    };

    struct CondensedSequenceGraph
    {
        std::vector<cds::ResidueType> types;
        std::vector<std::string> derivatives;
        graph::Graph graph;
    };

    SequenceGraph toSequenceGraph(const std::vector<ParsedResidue*>& residues);
    SequenceGraph condensedSequenceGraph(const SequenceGraph& sequence);
} // namespace cdsCondensedSequence

#endif
