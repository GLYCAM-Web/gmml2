#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEPARSER_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/Graph/graphTypes.hpp"

#include <memory>
#include <string>
#include <vector>

namespace cdsCondensedSequence
{
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

    SequenceData parseSequence(std::string sequence);
} // namespace cdsCondensedSequence
#endif
