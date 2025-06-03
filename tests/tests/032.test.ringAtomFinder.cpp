#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/SubgraphMatching.hpp"
#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/TotalCycleDecomposition.hpp"

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/ConnectivityIdentifier.hpp"
#include "includes/MolecularModeling/TemplateGraph/LazyPrints/LazyPrinters.hpp"

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"

#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <memory>

using namespace glygraph;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    //  Requirement: Find all atom rings in a pdb file.
    //  Format is bash-ready. Order of the rings is irrelevant. Numbers trace the ring rather than jumping around.
    //  Example output format:
    //  Ring_IDs=("1" "2" "3") ## just increment - this is just a counter
    //  Ring_Atoms=(
    //            ["1"]="1 2 3 4 5 6"
    //            ["2"]="21 22 23 24 25 26"
    //            ["3"]="31 32 33 34 35 36"
    //            )
    pdb::PdbFile pdbFile(argv[1], {pdb::InputType::modelsAsMolecules, false});
    std::cout << "Bonding atoms by distance for assembly" << std::endl;
    size_t atomCount = pdbFile.data.indices.atomCount;
    pdb::bondAtomsByDistance(pdbFile.data, codeUtils::indexVector(atomCount));
    Graph<cds::Atom>* graph1 = new Graph<cds::Atom>(pdbFile.data.objects.atoms.at(0));

    id_connectivity::identifyConnectivity(*graph1);
    // connectivity checking my dude
    //    std::set<Edge<cds::Atom>*> unknownEdges;
    //    std::set<Edge<cds::Atom>*> leafEdges;
    //    std::set<Edge<cds::Atom>*> bridgeEdges;
    //    std::set<Edge<cds::Atom>*> cycleEdges;
    //
    //    std::set<Node<cds::Atom>*> unknownNodes;
    //    std::set<Node<cds::Atom>*> leafNodes;
    //    std::set<Node<cds::Atom>*> bridgeNodes;
    //    std::set<Node<cds::Atom>*> cycleNodes;
    //
    //    for (Node<cds::Atom>* currNode : graph1->getNodes())
    //    {
    //        for (Edge<cds::Atom>* currEdge : currNode->getEdges())
    //        {
    //            switch (currEdge->getConnectivityTypeIdentifier())
    //            {
    //                case ConnectivityType::UNKNOWN:
    //                    unknownEdges.insert(currEdge);
    //                    break;
    //                case ConnectivityType::BRIDGE:
    //                    bridgeEdges.insert(currEdge);
    //                    break;
    //                case ConnectivityType::LEAF:
    //                    leafEdges.insert(currEdge);
    //                    break;
    //                case ConnectivityType::INCYCLE:
    //                    cycleEdges.insert(currEdge);
    //                    break;
    //                default:
    //                    badBehavior(__LINE__, __func__, "couldn't get edge connection type");
    //            }
    //        }
    //        switch (currNode->getConnectivityTypeIdentifier())
    //        {
    //            case ConnectivityType::UNKNOWN:
    //                unknownNodes.insert(currNode);
    //                break;
    //            case ConnectivityType::BRIDGE:
    //                bridgeNodes.insert(currNode);
    //                break;
    //            case ConnectivityType::LEAF:
    //                leafNodes.insert(currNode);
    //                break;
    //            case ConnectivityType::INCYCLE:
    //                cycleNodes.insert(currNode);
    //                break;
    //            default:
    //                badBehavior(__LINE__, __func__, "couldn't get node connection type");
    //        }
    //    }

    std::function<std::vector<size_t>(const std::pair<std::unordered_set<glygraph::Node<cds::Atom>*>,
                                                      std::unordered_set<glygraph::Edge<cds::Atom>*>>&)>
        cycleAtomIds = [&](const std::pair<std::unordered_set<glygraph::Node<cds::Atom>*>,
                                           std::unordered_set<glygraph::Edge<cds::Atom>*>>& cycle)
    {
        const std::unordered_set<glygraph::Node<cds::Atom>*>& cycleAtoms = cycle.first;
        std::vector<size_t> result;
        result.reserve(cycleAtoms.size());
        for (auto& a : cycleAtoms)
        {
            result.push_back(codeUtils::indexOf(pdbFile.data.objects.atoms, a->getDerivedClass()));
        }
        return codeUtils::sorted(result);
    };

    graph::Graph traversableGraph = graph::identity(pdbFile.data.atomGraph);

    std::function<std::vector<size_t>(const std::vector<size_t>&)> atomIdsInCycleOrder =
        [&](const std::vector<size_t>& atomIds)
    {
        std::vector<size_t> result;
        result.reserve(atomIds.size());
        std::vector<bool> partOfCycle = codeUtils::indicesToBools(atomCount, atomIds);
        size_t current                = atomIds[0];
        result.push_back(current);
        std::vector<bool> traversed = codeUtils::indicesToBools(atomCount, {current});
        auto untraversedNeighbor    = [&](size_t atomId)
        {
            return partOfCycle[atomId] && !traversed[atomId];
        };
        while (result.size() < atomIds.size())
        {
            const std::vector<size_t>& neighbors = traversableGraph.nodes.nodeAdjacencies[current];
            auto it = std::find_if(neighbors.begin(), neighbors.end(), untraversedNeighbor);
            if (it == neighbors.end())
            {
                throw std::runtime_error(
                    "Error in ring traversal. This shouldn't happen, since we already know it's a ring");
            }
            current = *it;
            result.push_back(current);
            traversed[current] = true;
        }
        return result;
    };

    std::function<std::string(const uint&)> uintToString = [&](uint a)
    {
        return std::to_string(a);
    };

    std::vector<std::string> atomNumberStrings = codeUtils::vectorMap(uintToString, pdbFile.data.atoms.numbers);
    std::vector<std::vector<size_t>> g1Cycles =
        codeUtils::vectorMap(cycleAtomIds, cycle_decomp::totalCycleDetect(*graph1));
    std::vector<std::vector<size_t>> orderedCycles = codeUtils::vectorMap(atomIdsInCycleOrder, g1Cycles);
    std::string separator                          = "";
    std::cout << "Ring_IDs=(";
    for (int i = 0; i < g1Cycles.size(); ++i)
    {
        std::cout << separator << "\"" << i << "\"";
        separator = " ";
    }
    std::cout << ")\n";
    int prettyCounter = 0;
    std::cout << "Ring_Atoms=(\n";
    for (auto& atomIds : orderedCycles)
    {
        std::cout << "[\"" << prettyCounter << "\"]=\"";
        std::cout << codeUtils::join(" ", codeUtils::indicesToValues(atomNumberStrings, atomIds));
        std::cout << "\"\n";
        prettyCounter++;
    }
    std::cout << ")\n";

    return 0;
}
