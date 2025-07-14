#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceUtil.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/parsing.hpp"

#include <optional>
#include <vector>
#include <stdexcept>
#include <iostream>

namespace
{
    using cdsCondensedSequence::AbstractSequence;
    using cdsCondensedSequence::ChildPosition;
    using cdsCondensedSequence::DualPosition;
    using cdsCondensedSequence::ParentPosition;
    using cdsCondensedSequence::SequenceData;

    size_t addNode(SequenceData& result, const cdsCondensedSequence::ParsedResidueComponents& components)
    {
        size_t id = graph::addNode(result.graph);
        result.residues.fullString.push_back(components.fullString);
        result.residues.type.push_back(components.type);
        result.residues.name.push_back(components.name);
        result.residues.ringType.push_back(components.ringType);
        result.residues.configuration.push_back(components.configuration);
        result.residues.isomer.push_back(components.isomer);
        result.residues.preIsomerModifier.push_back(components.preIsomerModifier);
        result.residues.ringShape.push_back(components.ringShape);
        result.residues.modifier.push_back(components.modifier);
        result.residues.isInternal.push_back(false);
        result.residues.isDerivative.push_back(
            codeUtils::contains({cds::ResidueType::Deoxy, cds::ResidueType::Derivative}, components.type));
        return id;
    }

    void addEdge(SequenceData& result, size_t parent, size_t child, const std::string& name,
                 const cdsCondensedSequence::EdgePosition& position)
    {
        graph::addEdge(result.graph, {parent, child});
        result.edges.names.push_back(name);
        result.edges.position.push_back(position);
        result.residues.isInternal[parent] =
            result.residues.isInternal[parent] || (result.residues.type[child] != cds::ResidueType::Deoxy);
    }

    size_t instantiateAglycone(SequenceData& result, const cdsCondensedSequence::AglyconeNode& node)
    {
        return addNode(result, node.components);
    }

    size_t instantiateDerivative(SequenceData& result, size_t parent, const cdsCondensedSequence::DerivativeNode& node)
    {
        const std::string& linkage = node.components.linkage;
        std::optional<uint> pos    = codeUtils::parseUint(linkage.substr(0, 1)).value();
        size_t id                  = addNode(result, node.components);
        addEdge(result, parent, id, linkage, ParentPosition {pos.value()});
        return parent;
    }

    size_t instantiateDerivativeList(SequenceData& result, size_t parent, const AbstractSequence& data,
                                     const cdsCondensedSequence::DerivativeListNode& node)
    {
        for (size_t n : node.constituents)
        {
            instantiateNode(result, parent, data, n);
        }
        return parent;
    }

    size_t instantiateMonosaccharide(SequenceData& result, size_t parent, const AbstractSequence& data,
                                     const cdsCondensedSequence::MonosaccharideNode& node)
    {
        const std::string& linkage        = node.components.linkage;
        const std::string linkageName     = node.components.configuration + linkage;
        const cds::ResidueType parentType = result.residues.type[parent];
        bool noParent                     = parentType == cds::ResidueType::Undefined;
        bool isParentAglycone             = parentType == cds::ResidueType::Aglycone;
        bool isParentSugar                = parentType == cds::ResidueType::Sugar;
        size_t id                         = addNode(result, node.components);
        if (noParent)
        {}
        else if (isParentAglycone)
        {
            std::optional<uint> position = codeUtils::parseUint(linkage.substr(0, 1)).value();
            addEdge(result, parent, id, linkageName, ChildPosition {position.value()});
        }
        else if (isParentSugar)
        {
            std::optional<uint> firstPosition = codeUtils::parseUint(linkage.substr(0, 1)).value();
            std::optional<uint> secondPosition =
                (linkage.length() > 2) ? codeUtils::parseUint(linkage.substr(2, 1)).value() : std::optional<uint> {};
            if (!(secondPosition.has_value()))
            {
                throw std::runtime_error("Either linkage or parent residue must contain parent linkage position");
            }
            uint childPosition  = firstPosition.value();
            uint parentPosition = secondPosition.value();
            addEdge(result, parent, id, linkageName, DualPosition {childPosition, parentPosition});
        }
        else
        {
            throw std::runtime_error("Parent of sugar must be aglycone or other sugar");
        }
        instantiateNode(result, id, data, node.derivativeList);
        return id;
    }

    size_t instantiateBranch(SequenceData& result, size_t parent, const AbstractSequence& data,
                             const cdsCondensedSequence::BranchNode& node)
    {
        instantiateNode(result, parent, data, node.chain);
        return parent;
    }

    size_t instantiateChain(SequenceData& result, size_t parent, const AbstractSequence& data,
                            const cdsCondensedSequence::ChainNode& node)
    {
        for (size_t n : node.constituents)
        {
            parent = instantiateNode(result, parent, data, n);
        }
        return parent;
    }
} // namespace

size_t cdsCondensedSequence::instantiateNode(SequenceData& result, size_t parent, const AbstractSequence& data,
                                             size_t id)
{
    const SequenceNode& node = data.nodes[id];
    if (std::holds_alternative<MonosaccharideNode>(node))
    {
        return instantiateMonosaccharide(result, parent, data, std::get<MonosaccharideNode>(node));
    }
    else if (std::holds_alternative<AglyconeNode>(node))
    {
        return instantiateAglycone(result, std::get<AglyconeNode>(node));
    }
    else if (std::holds_alternative<DerivativeNode>(node))
    {
        return instantiateDerivative(result, parent, std::get<DerivativeNode>(node));
    }
    else if (std::holds_alternative<DerivativeListNode>(node))
    {
        return instantiateDerivativeList(result, parent, data, std::get<DerivativeListNode>(node));
    }
    else if (std::holds_alternative<ChainNode>(node))
    {
        return instantiateChain(result, parent, data, std::get<ChainNode>(node));
    }
    else if (std::holds_alternative<BranchNode>(node))
    {
        return instantiateBranch(result, parent, data, std::get<BranchNode>(node));
    }
    else
    {
        throw std::runtime_error("Unhandled node type");
    }
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::instantiate(const AbstractSequence& data)
{
    SequenceData result;
    // use an empty node as parent when adding root residue
    addNode(result, {"", cds::ResidueType::Undefined, "", "", "", "", "", "", "", ""});
    result.graph.nodeAlive[0] = false;
    instantiateNode(result, 0, data, data.root);
    return rearrange(result, codeUtils::boolsToIndices(result.graph.nodeAlive), result.graph.edges);
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::rearrange(const SequenceData& sequence,
                                                                   const std::vector<size_t>& residueOrder,
                                                                   const std::vector<size_t>& edgeOrder)
{
    size_t residueCount = nodeCount(sequence.graph);
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

    for (size_t n : edgeOrder)
    {
        const std::array<size_t, 2>& edge = sequence.graph.edgeNodes[n];
        graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
    }

    ResidueData residues = {codeUtils::indicesToValues(sequence.residues.fullString, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.type, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.name, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringType, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.configuration, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isomer, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.preIsomerModifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.ringShape, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.modifier, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isInternal, residueOrder),
                            codeUtils::indicesToValues(sequence.residues.isDerivative, residueOrder)};

    EdgeData edges = {codeUtils::indicesToValues(sequence.edges.names, edgeOrder),
                      codeUtils::indicesToValues(sequence.edges.position, edgeOrder)};

    return SequenceData {resultGraph, residues, edges};
}

std::vector<size_t> cdsCondensedSequence::edgesSortedByLink(const SequenceData& sequence,
                                                            const std::vector<size_t>& edgeIds)
{
    std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
    {
        if (sequence.graph.edgeNodes[n][0] == sequence.graph.edgeNodes[k][0])
        {
            return mainLinkage(sequence, sequence.graph.edgeNodes[n][1]) >
                   mainLinkage(sequence, sequence.graph.edgeNodes[k][1]);
        }
        else
        {
            return sequence.graph.edgeNodes[n][0] < sequence.graph.edgeNodes[k][0];
        }
    };

    return codeUtils::sortedBy(compare, edgeIds);
}

cdsCondensedSequence::SequenceData cdsCondensedSequence::reordered(const SequenceData& sequence)
{
    size_t residueCount           = nodeCount(sequence.graph);
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

    return rearrange(sequence, residueOrder, edgeOrder);
}
