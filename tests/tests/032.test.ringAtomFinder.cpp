#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/SubgraphMatching.hpp"
#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/TotalCycleDecomposition.hpp"

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"

#include "includes/MolecularModeling/TemplateGraph/LazyPrints/LazyPrinters.hpp"
#include "includes/CodeUtils/casting.hpp"

#include "includes/MolecularModeling/TemplateGraph/Algorithms/include/ConnectivityIdentifier.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"

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
    pdb::bondAtomsByDistance(pdbFile.data, codeUtils::indexVector(pdbFile.data.indices.atomCount));
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

    // We had trouble casting the Node<cds::Atom>* to cds::Atom so came up with this instead.
    auto getAtomNumber = [&](Node<cds::Atom>* currDuder) -> unsigned int
    {
        size_t i = codeUtils::indexOf(pdbFile.data.objects.atoms, currDuder->getDerivedClass());
        return pdbFile.data.objects.atoms.at(i)->getNumber();
    };

    //    std::cout << "\nUnknown Nodes: ";
    //    for (Node<cds::Atom>* currDude : unknownNodes)
    //    {
    //    	std::cout << getAtomNumber(currDude) << ", ";
    //    }
    //    std::cout << "\nLeaf Nodes: ";
    //    for (Node<cds::Atom>* currDude : leafNodes)
    //    {
    //    	std::cout << getAtomNumber(currDude) << ", ";
    //    }
    //    std::cout << "\nBridge Nodes: ";
    //    for (Node<cds::Atom>* currDude : bridgeNodes)
    //    {
    //    	std::cout << getAtomNumber(currDude) << ", ";
    //    }
    //    std::cout << "\nCycle Nodes: ";
    //    for (Node<cds::Atom>* currDude : cycleNodes)
    //    {
    //    	std::cout << getAtomNumber(currDude) << ", ";
    //    }
    //    std::cout << "\n\n";

    std::vector<
        std::pair<std::unordered_set<glygraph::Node<cds::Atom>*>, std::unordered_set<glygraph::Edge<cds::Atom>*>>>
        g1Cycles          = cycle_decomp::totalCycleDetect(*graph1);
    std::string separator = "";
    std::cout << "Ring_IDs=(";
    for (int i = 0; i < g1Cycles.size(); ++i)
    {
        std::cout << separator << "\"" << i << "\"";
        separator = " ";
    }
    std::cout << ")\n";
    int prettyCounter = 0;
    std::cout << "Ring_Atoms=(\n";
    for (std::pair<std::unordered_set<glygraph::Node<cds::Atom>*>, std::unordered_set<glygraph::Edge<cds::Atom>*>>
             currCyclePair : g1Cycles)
    {
        std::string separator = "";
        std::cout << "[\"" << prettyCounter << "\"]=\"";
        for (Node<cds::Atom>* currAtom : currCyclePair.first)
        {
            std::cout << separator << getAtomNumber(currAtom);
            separator = " ";
        }
        std::cout << "\"\n";
        prettyCounter++;
    }
    std::cout << ")\n";

    return 0;
}
