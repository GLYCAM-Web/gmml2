#include "include/sequence/sequenceUtil.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace gmml
{
    namespace sequence
    {
        size_t parentEdge(const SequenceData& sequence, size_t nodeId)
        {
            auto childEquals = [&](const std::array<size_t, 2>& nodes) { return nodes[1] == nodeId; };
            auto it = std::find_if(sequence.graph.edges.nodes.begin(), sequence.graph.edges.nodes.end(), childEquals);
            return it - sequence.graph.edges.nodes.begin();
        }

        std::string edgeLinkage(const SequenceData& sequence, size_t edgeId)
        {
            const EdgePosition& position = sequence.edges.position[edgeId];
            if (std::holds_alternative<ParentPosition>(position))
            {
                return std::to_string(std::get<ParentPosition>(position).parentPosition);
            }
            else if (std::holds_alternative<ChildPosition>(position))
            {
                return std::to_string(std::get<ChildPosition>(position).childPosition) + "-";
            }
            else if (std::holds_alternative<DualPosition>(position))
            {
                const DualPosition& pos = std::get<DualPosition>(position);
                return std::to_string(pos.childPosition) + "-" + std::to_string(pos.parentPosition);
            }
            else
            {
                throw std::runtime_error("Unknown edge type");
            }
        }

        std::string mainLinkage(ResidueType type, const std::string& linkage)
        {
            switch (type)
            {
                case (ResidueType::Sugar):
                    return linkage.substr(linkage.size() - 1, 1);
                case (ResidueType::Derivative):
                    return linkage.substr(0, 1);
                case (ResidueType::Deoxy):
                    return linkage.substr(0, 1);
                default:
                    return std::string("0");
            }
        }

        std::string mainLinkage(const SequenceData& sequence, size_t residueId)
        {
            size_t edge = parentEdge(sequence, residueId);
            return mainLinkage(
                sequence.residues.type[residueId],
                (edge < edgeCount(sequence.graph)) ? edgeLinkage(sequence, edge) : "");
        }
    } // namespace sequence
} // namespace gmml
