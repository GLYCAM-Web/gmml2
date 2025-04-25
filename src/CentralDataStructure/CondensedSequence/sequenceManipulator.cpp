#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <sstream>
#include <ostream>
#include <functional>

using cdsCondensedSequence::ParsedResidue;

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    auto residueLink = [&](size_t n)
    {
        return cdsCondensedSequence::getLink(sequence.residues.type[n], sequence.residues.linkage[n]);
    };
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        return residueLink(sequence.graph.edgeNodes[n][1]) > residueLink(sequence.graph.edgeNodes[k][1]);
    };

    return codeUtils::sortedBy(compare, edgeIds);
}

std::vector<ParsedResidue*>
cdsCondensedSequence::parsedResiduesOrderedByConnectivity(std::vector<ParsedResidue*> residues)
{
    std::vector<ParsedResidue*> rawResidues;
    // Go via Graph so order decided by connectivity, depth first traversal:
    glygraph::Graph<cds::Residue> sequenceGraph(terminalResidue(residues));
    for (auto& node : sequenceGraph.getNodes())
    {
        rawResidues.push_back(codeUtils::erratic_cast<ParsedResidue*>(node->getDerivedClass()));
    }
    return rawResidues;
}

void cdsCondensedSequence::setIndexByConnectivity(std::vector<ParsedResidue*> residues)
{
    unsigned long long linkIndex    = 0; // Convention to start form 0 for linkages.
    unsigned long long residueIndex = 1; // Convention to start from 1 for residues.
    for (auto& residue : parsedResiduesOrderedByConnectivity(residues))
    {
        residue->setIndex(residueIndex);
        residue->setNumber(residueIndex); // ToDo temporary, switch to using number here. Keep index as a gmml internal
                                          // thing, never shown to user.
        ++residueIndex;
        for (auto& edge : residue->getInEdges())
        {
            edge->setIndex(linkIndex);
            ++linkIndex;
        }
    }
    return;
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::reordered(const SequenceData& sequence)
{
    size_t residueCount           = sequence.residues.name.size();
    std::vector<size_t> edgeOrder = edgesSortedByLink(sequence, codeUtils::indexVector(sequence.graph.edges));
    std::vector<std::array<size_t, 2>> reorderedEdges = codeUtils::indicesToValues(sequence.graph.edgeNodes, edgeOrder);

    size_t current                   = 0;
    std::vector<size_t> residueOrder = {current};
    residueOrder.reserve(residueCount);
    std::vector<bool> traversed = codeUtils::indicesToBools(residueCount, residueOrder);
    auto fromNode               = [&](size_t node)
    {
        std::vector<size_t> result;
        result.reserve(reorderedEdges.size());
        for (auto& edge : sequence.graph.edgeNodes)
        {
            if (edge[1] == node)
            {
                result.push_back(edge[0]);
            }
        }
        for (auto& edge : reorderedEdges)
        {
            if (edge[0] == node)
            {
                result.push_back(edge[1]);
            }
        }
        return codeUtils::reverse(result);
    };
    std::vector<size_t> nodesToTraverse;
    codeUtils::insertInto(nodesToTraverse, fromNode(current));
    while (!nodesToTraverse.empty())
    {
        current = nodesToTraverse.back();
        nodesToTraverse.pop_back();
        if (!traversed[current])
        {
            residueOrder.push_back(current);
            codeUtils::insertInto(nodesToTraverse, fromNode(current));
            traversed[current] = true;
        }
    }

    std::vector<size_t> missed = codeUtils::boolsToIndices(codeUtils::vectorNot(traversed));
    if (missed.size() > 0)
    {
        throw std::runtime_error("Error: sequence graph not fully connected");
    }

    std::vector<size_t> invertedResidueOrder(residueCount, -1);
    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        invertedResidueOrder[residueOrder[n]] = n;
    }

    graph::Database resultGraph;

    for (size_t n = 0; n < residueOrder.size(); n++)
    {
        graph::addNode(resultGraph);
    }

    std::vector<size_t> edgeIndexOrder = codeUtils::indexVector(sequence.graph.edges);

    for (size_t n : edgeIndexOrder)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    ResidueData residues = {codeUtils::indicesToValues(sequence.residues.fullString, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.type, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.name, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.linkage, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringType, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.configuration, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isomer, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringShape, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.modifier, residueOrder)};

    return SequenceData {resultGraph, residues, sequence.edges};
}

void cdsCondensedSequence::createParsedResidues(
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs, const SequenceData& sequence)
{
    size_t residueCount = sequence.residues.name.size();
    residuePtrs.reserve(residueCount);
    for (size_t n = 0; n < residueCount; n++)
    {
        residuePtrs.emplace_back(
            std::make_unique<cdsCondensedSequence::ParsedResidue>(cdsCondensedSequence::ParsedResidueComponents {
                sequence.residues.fullString[n],
                sequence.residues.type[n],
                sequence.residues.name[n],
                sequence.residues.linkage[n],
                sequence.residues.ringType[n],
                sequence.residues.configuration[n],
                sequence.residues.isomer[n],
                {sequence.residues.ringShape[n], sequence.residues.modifier[n]}
        }));
    }
    for (size_t n = 0; n < sequence.graph.edgeNodes.size(); n++)
    {
        auto& edge = sequence.graph.edgeNodes[n];
        residuePtrs[edge[1]].get()->addParent(sequence.edges.names[n], residuePtrs[edge[0]].get());
    }
}

void cdsCondensedSequence::sortResidueEdges(
    std::vector<std::unique_ptr<cdsCondensedSequence::ParsedResidue>>& residuePtrs)
{
    for (auto& residue : residuePtrs)
    {
        residue.get()->sortOutEdgesBySourceTObjectComparator();
    }
}

ParsedResidue* cdsCondensedSequence::terminalResidue(std::vector<ParsedResidue*> residues)
{
    return residues.front();
}
