#include "include/sequence/sequenceManipulation.hpp"

#include "include/graph/graphManipulation.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/sequence/sequenceUtil.hpp"
#include "include/util/containers.hpp"
#include "include/util/parsing.hpp"

#include <optional>
#include <stdexcept>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        namespace
        {
            size_t addNode(SequenceAndLinkageData& result, const ParsedResidueComponents& components)
            {
                ResidueData& residues = result.data.residues;
                size_t id = graph::addNode(result.data.graph);
                residues.fullString.push_back(components.fullString);
                residues.type.push_back(components.type);
                residues.name.push_back(components.name);
                residues.ringType.push_back(components.ringType);
                residues.configuration.push_back(components.configuration);
                residues.isomer.push_back(components.isomer);
                residues.preIsomerModifier.push_back(components.preIsomerModifier);
                residues.ringShape.push_back(components.ringShape);
                residues.modifier.push_back(components.modifier);
                residues.defaultHeadPosition.push_back(util::parseUint(components.defaultHeadPosition));
                residues.isInternal.push_back(false);
                residues.isDerivative.push_back(
                    util::contains({ResidueType::Deoxy, ResidueType::Derivative}, components.type));
                result.linkage.push_back(components.linkage);
                return id;
            }

            void addEdge(
                SequenceAndLinkageData& result,
                size_t parent,
                size_t child,
                const std::string& name,
                const EdgePosition& position)
            {
                ResidueData& residues = result.data.residues;
                EdgeData& edges = result.data.edges;
                graph::addEdge(result.data.graph, {parent, child});
                edges.names.push_back(name);
                edges.position.push_back(position);
                residues.isInternal[parent] =
                    residues.isInternal[parent] || (residues.type[child] != ResidueType::Deoxy);
            }

            void addMonosaccharideEdge(
                SequenceAndLinkageData& result,
                size_t parent,
                size_t child,
                const std::string& linkage,
                const std::string& configuration,
                const std::optional<uint>& defaultHeadPosition)
            {
                const ResidueType parentType = result.data.residues.type[parent];
                bool noParent = parentType == ResidueType::Undefined;
                bool isParentAglycone = parentType == ResidueType::Aglycone;
                bool isParentSugar = parentType == ResidueType::Sugar;
                if (noParent)
                {}
                else if (isParentAglycone)
                {
                    std::optional<uint> position = util::parseUint(linkage.substr(0, 1)).value();
                    addEdge(result, parent, child, configuration + linkage, ChildPosition {position.value()});
                }
                else if (isParentSugar)
                {
                    std::optional<uint> firstPosition = util::parseUint(linkage.substr(0, 1)).value();
                    std::optional<uint> secondPosition =
                        (linkage.length() > 2) ? util::parseUint(linkage.substr(2, 1)).value() : std::optional<uint> {};
                    if (!(secondPosition.has_value() || defaultHeadPosition.has_value()))
                    {
                        throw std::runtime_error(
                            "Either linkage or parent residue must contain parent linkage position");
                    }
                    uint childPosition = firstPosition.value();
                    uint parentPosition =
                        secondPosition.has_value() ? secondPosition.value() : defaultHeadPosition.value();
                    std::string linkageName =
                        configuration + std::to_string(childPosition) + "-" + std::to_string(parentPosition);
                    addEdge(result, parent, child, linkageName, DualPosition {childPosition, parentPosition});
                }
            }

            ChainState instantiateAglycone(SequenceAndLinkageData& result, const AglyconeNode& node)
            {
                size_t id = addNode(result, node.components);
                return {id, id};
            }

            ChainState instantiateDerivative(
                SequenceAndLinkageData& result, const ChainState& state, const DerivativeNode& node)
            {
                const std::string& linkage = node.components.linkage;
                std::optional<uint> pos = util::parseUint(linkage.substr(0, 1)).value();
                size_t id = addNode(result, node.components);
                size_t parent = state.head;
                addEdge(result, parent, id, linkage, ParentPosition {pos.value()});
                return state;
            }

            ChainState instantiateDerivativeList(
                SequenceAndLinkageData& result,
                const ChainState& state,
                const AbstractSequence& data,
                const DerivativeListNode& node)
            {
                for (size_t n : node.constituents)
                {
                    instantiateNode(result, state, data, n);
                }
                return state;
            }

            ChainState instantiateMonosaccharide(
                SequenceAndLinkageData& result,
                const ChainState& state,
                const AbstractSequence& data,
                const MonosaccharideNode& node)
            {
                size_t id = addNode(result, node.components);
                size_t parent = state.head;
                addMonosaccharideEdge(
                    result,
                    parent,
                    id,
                    node.components.linkage,
                    node.components.configuration,
                    result.data.residues.defaultHeadPosition[parent]);
                ChainState selfState {id, id};
                instantiateNode(result, selfState, data, node.derivativeList);
                return selfState;
            }

            ChainState instantiateChain(
                SequenceAndLinkageData& result,
                const ChainState& state,
                const AbstractSequence& data,
                const ChainNode& node)
            {
                if (node.constituents.empty())
                {
                    throw std::runtime_error("Chain cannot be empty");
                }
                ChainState first = instantiateNode(result, state, data, node.constituents[0]);
                ChainState current = first;
                for (size_t n = 1; n < node.constituents.size(); n++)
                {
                    current = instantiateNode(result, current, data, node.constituents[n]);
                }
                return {first.tail, current.head};
            }

            ChainState instantiateBranch(
                SequenceAndLinkageData& result,
                const ChainState& state,
                const AbstractSequence& data,
                const BranchNode& node)
            {
                ChainState empty {0, 0};
                ChainState current = instantiateNode(result, empty, data, node.chain);
                size_t parent = state.head;
                size_t tail = current.tail;
                std::optional<uint> position = util::parseUint(node.position);
                addMonosaccharideEdge(
                    result, parent, tail, result.linkage[tail], result.data.residues.configuration[tail], position);
                return state;
            }

            ChainState instantiateRepeat(
                SequenceAndLinkageData& result,
                const ChainState& state,
                const AbstractSequence& data,
                const RepeatNode& node)
            {
                ChainState first = instantiateNode(result, state, data, node.chain);
                ChainState current = first;
                for (size_t n = 1; n < node.repeats; n++)
                {
                    current = instantiateNode(result, current, data, node.chain);
                }
                return {first.tail, current.head};
            }
        } // namespace

        ChainState instantiateNode(
            SequenceAndLinkageData& result, const ChainState& state, const AbstractSequence& data, size_t id)
        {
            const SequenceNode& node = data.nodes[id];
            if (std::holds_alternative<MonosaccharideNode>(node))
            {
                return instantiateMonosaccharide(result, state, data, std::get<MonosaccharideNode>(node));
            }
            else if (std::holds_alternative<AglyconeNode>(node))
            {
                return instantiateAglycone(result, std::get<AglyconeNode>(node));
            }
            else if (std::holds_alternative<DerivativeNode>(node))
            {
                return instantiateDerivative(result, state, std::get<DerivativeNode>(node));
            }
            else if (std::holds_alternative<DerivativeListNode>(node))
            {
                return instantiateDerivativeList(result, state, data, std::get<DerivativeListNode>(node));
            }
            else if (std::holds_alternative<ChainNode>(node))
            {
                return instantiateChain(result, state, data, std::get<ChainNode>(node));
            }
            else if (std::holds_alternative<BranchNode>(node))
            {
                return instantiateBranch(result, state, data, std::get<BranchNode>(node));
            }
            else if (std::holds_alternative<RepeatNode>(node))
            {
                return instantiateRepeat(result, state, data, std::get<RepeatNode>(node));
            }
            else
            {
                throw std::runtime_error("Unhandled node type");
            }
        }

        SequenceData instantiate(const AbstractSequence& data)
        {
            SequenceAndLinkageData result;
            graph::Database& resultGraph = result.data.graph;
            // use an empty node as parent when adding root residue
            ParsedResidueComponents components;
            components.type = ResidueType::Undefined;
            addNode(result, components);
            resultGraph.nodes.alive[0] = false;
            instantiateNode(result, {0, 0}, data, data.root);
            if (!result.linkage[1].empty())
            {
                throw std::runtime_error(
                    "Sequence tail should not contain linkage in: " + result.data.residues.fullString[1]);
            }
            SequenceData rearr =
                rearrange(result.data, util::boolsToIndices(resultGraph.nodes.alive), resultGraph.edges.indices);
            return rearr;
        }

        SequenceData rearrange(
            const SequenceData& sequence, const std::vector<size_t>& residueOrder, const std::vector<size_t>& edgeOrder)
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
                const std::array<size_t, 2>& edge = sequence.graph.edges.nodes[n];
                graph::addEdge(resultGraph, {invertedResidueOrder[edge[0]], invertedResidueOrder[edge[1]]});
            }

            ResidueData residues = {
                util::indicesToValues(sequence.residues.fullString, residueOrder),
                util::indicesToValues(sequence.residues.type, residueOrder),
                util::indicesToValues(sequence.residues.name, residueOrder),
                util::indicesToValues(sequence.residues.ringType, residueOrder),
                util::indicesToValues(sequence.residues.configuration, residueOrder),
                util::indicesToValues(sequence.residues.isomer, residueOrder),
                util::indicesToValues(sequence.residues.preIsomerModifier, residueOrder),
                util::indicesToValues(sequence.residues.ringShape, residueOrder),
                util::indicesToValues(sequence.residues.modifier, residueOrder),
                util::indicesToValues(sequence.residues.defaultHeadPosition, residueOrder),
                util::indicesToValues(sequence.residues.isInternal, residueOrder),
                util::indicesToValues(sequence.residues.isDerivative, residueOrder)};

            EdgeData edges = {
                util::indicesToValues(sequence.edges.names, edgeOrder),
                util::indicesToValues(sequence.edges.position, edgeOrder)};

            return SequenceData {resultGraph, residues, edges};
        }

        std::vector<size_t> edgesSortedByLink(const SequenceData& sequence, const std::vector<size_t>& edgeIds)
        {
            std::function<bool(const size_t&, const size_t&)> compare = [&](const size_t& n, const size_t& k)
            {
                if (sequence.graph.edges.nodes[n][0] == sequence.graph.edges.nodes[k][0])
                {
                    return mainLinkage(sequence, sequence.graph.edges.nodes[n][1]) >
                           mainLinkage(sequence, sequence.graph.edges.nodes[k][1]);
                }
                else
                {
                    return sequence.graph.edges.nodes[n][0] < sequence.graph.edges.nodes[k][0];
                }
            };

            return util::sortedBy(compare, edgeIds);
        }

        SequenceData reordered(const SequenceData& sequence)
        {
            size_t residueCount = nodeCount(sequence.graph);
            std::vector<size_t> edgeOrder =
                edgesSortedByLink(sequence, util::indexVector(sequence.graph.edges.indices));
            std::vector<std::array<size_t, 2>> reorderedEdges =
                util::indicesToValues(sequence.graph.edges.nodes, edgeOrder);

            size_t current = 0;
            std::vector<size_t> residueOrder = {current};
            residueOrder.reserve(residueCount);
            std::vector<bool> traversed = util::indicesToBools(residueCount, residueOrder);
            auto fromNode = [&](size_t node)
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
                return util::reverse(result);
            };
            std::vector<size_t> nodesToTraverse;
            util::insertInto(nodesToTraverse, fromNode(current));
            while (!nodesToTraverse.empty())
            {
                current = nodesToTraverse.back();
                nodesToTraverse.pop_back();
                if (!traversed[current])
                {
                    residueOrder.push_back(current);
                    util::insertInto(nodesToTraverse, fromNode(current));
                    traversed[current] = true;
                }
            }

            std::vector<size_t> missed = util::boolsToIndices(util::vectorNot(traversed));
            if (missed.size() > 0)
            {
                throw std::runtime_error("Error: sequence graph not fully connected");
            }

            return rearrange(sequence, residueOrder, edgeOrder);
        }
    } // namespace sequence
} // namespace gmml
