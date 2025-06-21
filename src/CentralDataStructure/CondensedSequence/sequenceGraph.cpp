#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>
#include <array>
#include <vector>

std::vector<std::string> cdsCondensedSequence::sequenceDerivatives(const SequenceData& sequence)
{
    std::vector<std::string> derivatives(sequence.nodes.name.size(), "");

    for (size_t n = 0; n < sequence.graph.edgeNodes.size(); n++)
    {
        const std::array<size_t, 2>& adj = sequence.graph.edgeNodes[n];
        size_t parent                    = adj[0];
        size_t child                     = adj[1];
        if (sequence.nodes.type[child] == cds::ResidueType::Derivative)
        {
            std::string& str = derivatives[parent];
            if (!str.empty())
            {
                str += " ";
            }
            str += sequence.edges.names[n] + sequence.nodes.name[child];
        }
    }

    return derivatives;
}

std::vector<std::string> cdsCondensedSequence::sequenceMonosaccharideNames(const SequenceData& sequence)
{
    std::function<std::string(const size_t&)> monosaccharideName = [&](size_t n)
    {
        return sequence.nodes.isomer[n] + sequence.nodes.name[n] + sequence.nodes.modifier[n] +
               sequence.nodes.ringShape[n];
    };
    return codeUtils::vectorMap(monosaccharideName, codeUtils::indexVector(sequence.nodes.name.size()));
}

graph::Graph cdsCondensedSequence::condensedSequenceGraph(const SequenceData& sequence)
{
    std::function<bool(const cds::ResidueType&)> keepNode = [](const cds::ResidueType& type)
    {
        return type != cds::ResidueType::Derivative;
    };
    graph::Database db = sequence.graph;
    db.nodeAlive       = codeUtils::vectorMap(keepNode, sequence.nodes.type);
    return graph::identity(db);
}
