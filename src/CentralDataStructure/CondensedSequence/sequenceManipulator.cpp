#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <sstream>
#include <sys/stat.h> // for checking if file exists
#include <fstream>    // writing outputDotFile

using cdsCondensedSequence::ParsedResidue;

namespace
{
    bool file_exists(const char* filename)
    {
        struct stat buffer;
        return (stat(filename, &buffer) == 0);
    }

    std::string getGraphVizLineForResidue(ParsedResidue& residue, cdsCondensedSequence::GraphVizDotConfig& configs)
    {
        std::stringstream logss;
        logss << "Getting GraphVizLine for " << residue.GetName() << "\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
        logss.str(std::string());
        logss.clear(); // Must do both of these to clear the stream;
        std::stringstream ss;
        ss << residue.getIndex() << " [";
        // Aglycone
        if (residue.GetType() == cds::ResidueType::Aglycone)
        {
            ss << "shape=box label=\"" << residue.GetMonosaccharideName() << "\"]";
            return ss.str();
        }
        // Sugar
        std::string imageFile = configs.svg_directory_path_ + residue.GetMonosaccharideName() + ".svg";
        logss << "Searching for image: " << imageFile << "\n";
        if (file_exists(imageFile.c_str()))
        {
            logss << "FOUND IT\n";
            std::string label = (residue.GetRingType() == "f") ? "f" : "";
            ss << "label=\"" << label << "\" height=\"0.7\" image=\"" << imageFile << "\"];\n";
        }
        else
        {
            logss << "Not image available, using circle\n";
            ss << "shape=circle height=\"0.7\" label=\"" << residue.GetMonosaccharideName() << "\"];\n";
        }
        // Derivatives
        std::string derivativeStr = "";
        for (auto& childLink : residue.GetChildren())
        {
            if (childLink->GetType() == cds::ResidueType::Derivative)
            {
                derivativeStr += childLink->GetLinkageName() + childLink->GetName() + " ";
            }
        }
        if (!derivativeStr.empty())
        {
            ss << "\n"
               << "b" << residue.getIndex();
            ss << "[ shape=\"plaintext\",fontsize=\"12\",forcelabels=\"true\"; height = \"0.3\"; labelloc = b;  "
                  "label=\"";
            ss << derivativeStr << "\"];\n";
            ss << "{ rank=\"same\"; b" << residue.getIndex() << " " << residue.getIndex() << "};\n";
            ss << "{nodesep=\"0.2\";b" << residue.getIndex() << ";" << residue.getIndex() << "};\n";
            ss << "b" << residue.getIndex() << "--" << residue.getIndex() << " [style=invis];\n";
        }
        // Linkage
        for (auto& parent : residue.GetParents())
        { // There is either 1 or 0 parents, this covers both cases.
            ss << residue.getIndex() << "--" << parent->getIndex() << "[label=\"";
            if (configs.show_config_labels_)
            {
                ss << residue.GetConfiguration();
            }
            if (configs.show_position_labels_)
            {
                ss << residue.GetLinkage();
            }
            ss << "\"];\n";
            if (configs.show_edge_labels_)
            {
                for (auto& linkage : residue.getOutEdges())
                {
                    ss << residue.getIndex() << "--" << parent->getIndex();
                    ss << "[taillabel=< <B>" << linkage->getIndex() << "</B>>, ";
                    ss << "labelfontsize = 14, labeldistance = 2.0, labelangle = -35";
                    ss << "];\n";
                }
            }
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
        return ss.str();
    }

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
    setIndexByConnectivity(residues);
    std::stringstream ss;
    ss << "graph G {graph [splines=false forcelabels=true  dpi=" << configs.dpi_ << "];\n";
    ss << "node [ shape=\"none\" fontname=DejaVuSans labelfontsize=12 forcelabels=\"true\";\n";
    ss << "label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];\n";
    ss << "edge [labelfontsize=12 fontname=DejaVuSans labeldistance=1.2 labelangle = 320.0];\n";
    ss << "rankdir=LR nodesep=\"0.05\" ranksep=\"0.8\";\n";
    for (auto& residue : parsedResiduesOrderedByConnectivity(residues))
    {
        if (residue->GetType() != cds::ResidueType::Derivative)
        {
            ss << getGraphVizLineForResidue(*residue, configs) << "\n";
        }
    }
    ss << "}\n";
    // Open and overwrite.
    std::ofstream outputDotFile(configs.file_name_, std::ios::trunc);
    outputDotFile << ss.str();
    outputDotFile.close();
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return ss.str();
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
