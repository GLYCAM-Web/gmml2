#include "include/sequence/sequenceGraph.hpp"

#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/residueTypes.hpp"
#include "include/util/containers.hpp"

#include <array>
#include <functional>
#include <vector>

namespace gmml
{
    namespace sequence
    {
        std::vector<std::string> sequenceDerivatives(const SequenceData& sequence)
        {
            std::vector<std::string> derivatives(nodeCount(sequence.graph), "");

            for (size_t n = 0; n < edgeCount(sequence.graph); n++)
            {
                const std::array<size_t, 2>& adj = sequence.graph.edgeNodes[n];
                size_t parent = adj[0];
                size_t child = adj[1];
                if (sequence.residues.type[child] == ResidueType::Derivative)
                {
                    std::string& str = derivatives[parent];
                    if (!str.empty())
                    {
                        str += " ";
                    }
                    str += sequence.edges.names[n] + sequence.residues.name[child];
                }
            }

            return derivatives;
        }

        std::vector<std::string> sequenceMonosaccharideNames(const SequenceData& sequence)
        {
            std::function<std::string(const size_t&)> monosaccharideName = [&](size_t n)
            {
                return sequence.residues.isomer[n] + sequence.residues.name[n] + sequence.residues.modifier[n] +
                       sequence.residues.ringShape[n];
            };
            return util::vectorMap(monosaccharideName, util::indexVector(nodeCount(sequence.graph)));
        }

        graph::Graph condensedSequenceGraph(const SequenceData& sequence)
        {
            std::function<bool(const ResidueType&)> keepNode = [](const ResidueType& type)
            { return type != ResidueType::Derivative; };
            graph::Database db = sequence.graph;
            db.nodeAlive = util::vectorMap(keepNode, sequence.residues.type);
            return graph::identity(db);
        }
    } // namespace sequence
} // namespace gmml
