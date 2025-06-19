#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <memory>
#include <string>
#include <vector>
#include <variant>

namespace cdsCondensedSequence
{
    struct ResidueNode
    {};

    struct BranchNode
    {
        std::string linkage;
        size_t head;
    };

    struct ProbabilityNode
    {
        double probability;
        std::vector<size_t> heads;
        std::vector<size_t> tails;
        std::vector<double> weights;
    };

    typedef std::variant<ResidueNode, BranchNode, ProbabilityNode> NodeType;

    struct NodeData
    {
        std::vector<std::string> fullString;
        std::vector<cds::ResidueType> type;
        std::vector<std::string> name;
        std::vector<std::string> linkage;
        std::vector<std::string> chainPosition;
        std::vector<std::string> ringType;
        std::vector<std::string> configuration;
        std::vector<std::string> isomer;
        std::vector<std::string> preIsomerModifier;
        std::vector<std::string> ringShape;
        std::vector<std::string> modifier;
        std::vector<bool> isInternal;
        std::vector<bool> isDerivative;
        std::vector<NodeType> nodeType;
    };

    struct EdgeData
    {
        std::vector<std::string> names;
    };

    struct SequenceData
    {
        graph::Database graph;
        NodeData nodes;
        EdgeData edges;
    };

    SequenceData parseSequence(std::string sequence);
} // namespace cdsCondensedSequence
#endif
