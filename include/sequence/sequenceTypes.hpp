#ifndef INCLUDE_SEQUENCE_SEQUENCETYPES_HPP
#define INCLUDE_SEQUENCE_SEQUENCETYPES_HPP

#include "include/CentralDataStructure/residueTypes.hpp"
#include "include/graph/graphTypes.hpp"

#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        struct ParsedResidueComponents
        {
            std::string fullString;
            ResidueType type;
            std::string name;
            std::string linkage;
            std::string defaultHeadPosition;
            std::string ringType;
            std::string configuration;
            std::string isomer;
            std::string preIsomerModifier;
            std::string ringShape;
            std::string modifier;
        };

        struct MonosaccharideNode
        {
            ParsedResidueComponents components;
            size_t derivativeList;
        };

        struct AglyconeNode
        {
            ParsedResidueComponents components;
        };

        struct DerivativeNode
        {
            ParsedResidueComponents components;
        };

        struct DerivativeListNode
        {
            std::vector<size_t> constituents;
        };

        struct ChainNode
        {
            std::vector<size_t> constituents;
        };

        struct BranchNode
        {
            std::string position;
            size_t chain;
        };

        struct RepeatNode
        {
            uint repeats;
            size_t chain;
        };

        typedef std::variant<
            MonosaccharideNode,
            AglyconeNode,
            DerivativeNode,
            DerivativeListNode,
            ChainNode,
            BranchNode,
            RepeatNode>
            SequenceNode;

        struct AbstractSequence
        {
            std::vector<std::string> fullString;
            std::vector<SequenceNode> nodes;
            size_t root;
        };

        struct ResidueData
        {
            std::vector<std::string> fullString;
            std::vector<ResidueType> type;
            std::vector<std::string> name;
            std::vector<std::string> ringType;
            std::vector<std::string> configuration;
            std::vector<std::string> isomer;
            std::vector<std::string> preIsomerModifier;
            std::vector<std::string> ringShape;
            std::vector<std::string> modifier;
            std::vector<std::optional<uint>> defaultHeadPosition;
            std::vector<bool> isInternal;
            std::vector<bool> isDerivative;
        };

        struct ChildPosition
        {
            uint childPosition;
        };

        struct ParentPosition
        {
            uint parentPosition;
        };

        struct DualPosition
        {
            uint childPosition;
            uint parentPosition;
        };

        typedef std::variant<ChildPosition, ParentPosition, DualPosition> EdgePosition;

        struct EdgeData
        {
            std::vector<std::string> names;
            std::vector<EdgePosition> position;
        };

        struct SequenceData
        {

            graph::Database graph;
            ResidueData residues;
            EdgeData edges;
        };
    } // namespace sequence
} // namespace gmml

#endif
