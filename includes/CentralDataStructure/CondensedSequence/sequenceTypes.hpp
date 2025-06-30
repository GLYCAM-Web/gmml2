#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCETYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCETYPES_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/residueTypes.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <string>
#include <vector>

namespace cdsCondensedSequence
{
    struct ResidueNode
    {
        ParsedResidueComponents components;
    };

    struct AbstractSequence
    {
        graph::Database graph;
        std::vector<ResidueNode> nodes;
    };

    struct ResidueData
    {
        std::vector<std::string> fullString;
        std::vector<cds::ResidueType> type;
        std::vector<std::string> name;
        std::vector<std::string> linkage;
        std::vector<std::string> ringType;
        std::vector<std::string> configuration;
        std::vector<std::string> isomer;
        std::vector<std::string> preIsomerModifier;
        std::vector<std::string> ringShape;
        std::vector<std::string> modifier;
        std::vector<bool> isInternal;
        std::vector<bool> isDerivative;
    };

    struct EdgeData
    {
        std::vector<std::string> names;
    };

    struct SequenceData
    {
        graph::Database graph;
        ResidueData residues;
        EdgeData edges;
    };
} // namespace cdsCondensedSequence
#endif
