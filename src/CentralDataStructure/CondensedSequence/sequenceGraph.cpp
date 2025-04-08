#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <functional>
#include <array>
#include <vector>

cdsCondensedSequence::SequenceGraph cdsCondensedSequence::toSequenceGraph(const std::vector<ParsedResidue*>& residues)
{
    size_t residueCount = residues.size();
    std::vector<uint> initialIndices;
    initialIndices.reserve(residueCount);
    SequenceGraph result;
    result.types.reserve(residueCount);
    result.names.reserve(residueCount);
    result.monosaccharideNames.reserve(residueCount);
    result.ringTypes.reserve(residueCount);
    result.linkageNames.reserve(residueCount);
    result.configurations.reserve(residueCount);
    result.linkages.reserve(residueCount);
    result.derivatives = std::vector<std::string>(residueCount, "");
    graph::Database db;
    for (uint n = 0; n < residueCount; n++)
    {
        ParsedResidue* residue = residues[n];
        graph::addNode(db);
        initialIndices.push_back(residue->getIndex());
        residue->setIndex(n);
        result.types.push_back(residue->GetType());
        result.names.push_back(residue->GetName());
        std::string monosaccharideName =
            residue->GetIsomer() + residue->GetResidueName() + residue->GetResidueModifier() + residue->GetRingShape();
        result.monosaccharideNames.push_back(monosaccharideName);
        result.ringTypes.push_back(residue->GetRingType());
        result.linkageNames.push_back(residue->GetLinkageName());
        result.configurations.push_back(residue->GetConfiguration());
        result.linkages.push_back(residue->GetLinkage());
    }
    for (size_t n = 0; n < residueCount; n++)
    {
        for (auto& neighbor : residues[n]->getChildren())
        {
            uint index = neighbor->getIndex();
            addEdge(db, {n, index});
        }
    }
    for (size_t n = 0; n < residueCount; n++)
    {
        residues[n]->setIndex(initialIndices[n]);
    }
    result.graph = graph::identity(db);
    return result;
}

cdsCondensedSequence::SequenceGraph cdsCondensedSequence::condensedSequenceGraph(const SequenceGraph& sequence)
{
    std::function<bool(const cds::ResidueType&)> keepNode = [](const cds::ResidueType& type)
    {
        return type != cds::ResidueType::Derivative;
    };
    std::vector<std::string> derivatives = sequence.derivatives;
    for (auto& adj : sequence.graph.edges.nodeAdjacencies)
    {
        size_t parent = adj[0];
        size_t child  = adj[1];
        if (sequence.types[child] == cds::ResidueType::Derivative)
        {
            std::string& str = derivatives[parent];
            if (!str.empty())
            {
                str += " ";
            }
            str += sequence.linkageNames[child] + sequence.names[child];
        }
    }
    graph::Database db   = sequence.graph.source;
    db.nodeAlive         = codeUtils::vectorMap(keepNode, sequence.types);
    SequenceGraph result = sequence;
    result.derivatives   = derivatives;
    result.graph         = graph::identity(db);
    return result;
}
