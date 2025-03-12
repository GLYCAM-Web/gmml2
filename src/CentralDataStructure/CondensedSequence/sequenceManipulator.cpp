#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceGraph.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <sstream>
#include <ostream>

using cdsCondensedSequence::ParsedResidue;

namespace
{
    void recurvePrintIupac(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output)
    {
        auto neighbors                  = currentResidue->GetChildren();
        size_t numberOfNeighbors        = neighbors.size();
        // Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
        std::string outputResidueString = currentResidue->GetIupacName();
        std::vector<std::string> derivatives;
        for (auto& neighbor : neighbors)
        {
            if (neighbor->GetType() == cds::ResidueType::Derivative || neighbor->GetType() == cds::ResidueType::Deoxy)
            {
                --numberOfNeighbors;
                derivatives.push_back(neighbor->GetLinkageName() + neighbor->GetName());
                // derivatives.push_back(","); not for iupac it looks like Glc2Me3Ac
            }
        }
        if (!derivatives.empty())
        {
            // derivatives.pop_back();                               // Remove the last ","
            std::reverse(derivatives.begin(), derivatives.end()); // order should be 2S,6S, not 6S,2S.
            // outputResidueString += "[";
            for (auto& derivative : derivatives)
            {
                outputResidueString += derivative;
            }
            // outputResidueString += "]";
        }
        // Output
        if (currentResidue->GetType() != cds::ResidueType::Aglycone)
        { // needs () around the linkageName, but not the aglycone
            outputResidueString += "(";
        }
        outputResidueString += currentResidue->GetLinkageName();
        if (currentResidue->GetType() != cds::ResidueType::Aglycone)
        { // Reducing/rightmost residue has no parent.
            if (currentResidue->getParent() != nullptr &&
                currentResidue->getParent()->GetType() != cds::ResidueType::Aglycone)
            { // IUPAC leaves the terminal open. So you get Glc(b1-
                outputResidueString += ")";
            }
        }
        output.push_back(outputResidueString);
        // End of a branch check
        if (numberOfNeighbors == 0 && branchStackSize > 0)
        {
            output.push_back("[");
            --branchStackSize;
        }
        size_t loopCount = 0;
        for (auto& neighbor : neighbors)
        {
            if (neighbor->GetType() != cds::ResidueType::Derivative && neighbor->GetType() != cds::ResidueType::Deoxy)
            {
                ++loopCount;
                if (loopCount < numberOfNeighbors)
                {
                    output.push_back("]");
                    ++branchStackSize;
                }
                recurvePrintIupac(neighbor, branchStackSize, output);
            }
        }
        return;
    }

    void recurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output,
                      const bool withLabels)
    {
        auto neighbors                  = currentResidue->GetChildren();
        size_t numberOfNeighbors        = neighbors.size();
        // Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
        std::string outputResidueString = currentResidue->GetName(withLabels);
        std::vector<std::string> derivatives;
        for (auto& neighbor : neighbors)
        {
            if (neighbor->GetType() == cds::ResidueType::Derivative || neighbor->GetType() == cds::ResidueType::Deoxy)
            {
                --numberOfNeighbors;
                derivatives.push_back(neighbor->GetLinkageName(withLabels) + neighbor->GetName(withLabels));
                derivatives.push_back(",");
            }
        }
        if (!derivatives.empty())
        {
            derivatives.pop_back();                               // Remove the last ","
            std::reverse(derivatives.begin(), derivatives.end()); // order should be 2S,6S, not 6S,2S.
            outputResidueString += "[";
            for (auto& derivative : derivatives)
            {
                outputResidueString += derivative;
            }
            outputResidueString += "]";
        }
        outputResidueString += currentResidue->GetLinkageName(withLabels);
        output.push_back(outputResidueString);
        // End of a branch check
        if (numberOfNeighbors == 0 && branchStackSize > 0)
        {
            output.push_back("[");
            --branchStackSize;
        }
        size_t loopCount = 0;
        for (auto& neighbor : neighbors)
        {
            if (neighbor->GetType() != cds::ResidueType::Derivative && neighbor->GetType() != cds::ResidueType::Deoxy)
            {
                ++loopCount;
                if (loopCount < numberOfNeighbors)
                {
                    output.push_back("]");
                    ++branchStackSize;
                }
                recurvePrint(neighbor, branchStackSize, output, withLabels);
            }
        }
        return;
    }
} // namespace

std::string cdsCondensedSequence::printGraphViz(GraphVizDotConfig& configs, std::vector<ParsedResidue*> residues)
{
    SequenceGraph sequence = condensedSequenceGraph(toSequenceGraph(residues));

    std::vector<GraphVizResidueNode> nodes;
    nodes.reserve(nodeCount(sequence.graph));
    std::vector<GraphVizLinkage> linkages;
    linkages.reserve(edgeCount(sequence.graph));

    for (size_t n = 0; n < nodeCount(sequence.graph); n++)
    {
        size_t index                    = sourceNodeIndex(sequence.graph, n);
        std::string& monosaccharideName = sequence.monosaccharideNames[index];
        GraphVizImage image = findImage(configs, monosaccharideName, (sequence.ringTypes[index] == "f") ? "f" : "");
        nodes.push_back({n, image, monosaccharideName, sequence.derivatives[index]});
    }

    for (size_t n = 0; n < edgeCount(sequence.graph); n++)
    {
        std::array<size_t, 2>& adj = sequence.graph.edges.nodeAdjacencies[n];
        size_t parent              = adj[0];
        size_t child               = adj[1];
        size_t childIndex          = sourceNodeIndex(sequence.graph, child);
        std::string label          = (configs.show_config_labels_ ? sequence.configurations[childIndex] : "") +
                            (configs.show_position_labels_ ? sequence.linkages[childIndex] : "");
        linkages.push_back({
            {child, parent},
            label
        });
    }

    std::stringstream ss;
    ss << "graph G {graph [splines=false dpi=" << configs.dpi_ << " outputorder=\"edgesfirst\"];\n";
    ss << "node [shape=\"none\" fontname=DejaVuSans labelfontsize=12 ";
    ss << "label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];\n";
    ss << "edge [labelfontsize=12 fontname=DejaVuSans labeldistance=1.2 labelangle=320.0];\n";
    ss << "rankdir=LR nodesep=\"0.05\" ranksep=\"0.8\";\n";
    for (size_t n = 0; n < nodes.size(); n++)
    {
        size_t nodeIndex = sourceNodeIndex(sequence.graph, n);
        bool isAglycone  = sequence.types[nodeIndex] == cds::ResidueType::Aglycone;
        ss << (isAglycone ? graphVizAglyconeNode(nodes[n]) : graphVizSugarNode(nodes[n]));
    }
    for (auto& linkage : linkages)
    {
        ss << graphVizLinkageLine(nodes, linkage);
    }
    ss << "}\n";
    std::string str = ss.str();
    // Open and overwrite.
    codeUtils::writeToFile(configs.file_name_,
                           [&str](std::ostream& stream)
                           {
                               stream << str;
                           });
    return str;
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

void cdsCondensedSequence::labelSequence(std::vector<ParsedResidue*> residues)
{
    setIndexByConnectivity(residues);
    std::stringstream ss;
    for (auto& residue : parsedResiduesOrderedByConnectivity(residues))
    {
        ss << residue->GetName() << "&Label=residue-" << residue->getIndex() << ";";
        residue->addLabel(ss.str());
        ss.str(std::string());
        ss.clear(); // Must do both of these to clear the stream
        for (auto& edge : residue->getOutEdges())
        {
            ss << edge->getLabel() << "&Label=link-" << edge->getIndex() << ";";
            edge->addLabel(ss.str());
            ss.str(std::string());
            ss.clear(); // Must do both of these to clear the stream
        }
    }
    return;
}

std::string cdsCondensedSequence::printSequence(std::vector<ParsedResidue*> residues, bool withLabels,
                                                bool iupacCondensed)
{
    std::vector<std::string> output;
    int branchStackSize = 0;
    if (iupacCondensed)
    {
        recurvePrintIupac(terminalResidue(residues), branchStackSize, output);
    }
    else
    {
        recurvePrint(terminalResidue(residues), branchStackSize, output, withLabels);
    }
    std::reverse(output.begin(), output.end()); // Reverse order, as it starts from terminal.
    std::stringstream ss;
    for (auto& label : output)
    {
        ss << label;
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
}

std::string cdsCondensedSequence::reorderSequence(std::vector<std::unique_ptr<ParsedResidue>>& residues)
{ // Just doing the default by ascending link number for now.
    for (auto& residue : residues)
    {
        residue.get()->sortOutEdgesBySourceTObjectComparator();
    }
    glygraph::Graph<cds::Residue> sequenceGraph(residues.front().get());
    size_t newPosition = 0;
    for (auto& node : sequenceGraph.getNodes())
    {
        for (size_t oldPosition = 0; oldPosition < residues.size(); oldPosition++)
        {
            if (residues[oldPosition].get() == node->getDerivedClass())
            {
                if (oldPosition != newPosition)
                {
                    std::swap(residues[oldPosition], residues[newPosition]);
                }
                break;
            }
        }
        newPosition++;
    }
    return printSequence(codeUtils::pointerToUniqueVector(residues), false);
}

ParsedResidue* cdsCondensedSequence::terminalResidue(std::vector<ParsedResidue*> residues)
{
    return residues.front();
}
