#include "includes/CentralDataStructure/CondensedSequence/graphViz.hpp"
#include "includes/CodeUtils/files.hpp"

#include <string>
#include <array>
#include <vector>

cdsCondensedSequence::GraphVizImage
cdsCondensedSequence::findImage(const cdsCondensedSequence::GraphVizDotConfig& configs, const std::string& name,
                                const std::string& label)
{
    std::string imageFile = configs.svg_directory_path_ + name + ".svg";
    if (codeUtils::doesFileExist((configs.base_path_ + imageFile).c_str()))
    {
        return GraphVizImage {true, imageFile, label};
    }
    return GraphVizImage {false, "", ""};
}

std::string cdsCondensedSequence::graphVizAglyconeNode(const GraphVizResidueNode& node)
{
    std::stringstream ss;
    ss << node.index << " [shape=box label=\"" << node.label << "\"]\n";
    return ss.str();
}

std::string cdsCondensedSequence::graphVizSugarNode(const GraphVizResidueNode& node)
{
    std::stringstream ss;
    ss << node.index;
    if (node.image.found)
    {
        ss << " [label=\"" << node.image.label << "\" height=\"0.7\" image=\"" << node.image.path << "\"];\n";
    }
    else
    {
        ss << " [shape=circle height=\"0.7\" label=\"" << node.label << "\"];\n";
    }
    // Derivatives
    if (!node.floatingLabel.empty())
    {
        ss << "b" << node.index;
        ss << " [shape=\"plaintext\" fontsize=\"12\" height=\"0.3\" labelloc=b label=\"";
        ss << node.floatingLabel << "\"];\n";
        ss << "{rank=\"same\" b" << node.index << " " << node.index << "};\n";
        ss << "{nodesep=\"0.2\" b" << node.index << " " << node.index << "};\n";
        ss << "b" << node.index << "--" << node.index << " [style=invis];\n";
    }
    return ss.str();
}

std::string cdsCondensedSequence::graphVizLinkageLine(const std::vector<GraphVizResidueNode>& nodes,
                                                      const GraphVizLinkage& linkage)
{
    std::array<std::string, 2> clip;
    for (size_t n = 0; n < 2; n++)
    {
        clip[n] = nodes[linkage.nodeIndices[n]].image.found ? "false" : "true";
    }
    std::stringstream ss;
    ss << nodes[linkage.nodeIndices[0]].index << "--" << nodes[linkage.nodeIndices[1]].index << " [label=\""
       << linkage.label << "\" headclip=" << clip[1] << " tailclip=" << clip[0] << "];\n";
    return ss.str();
}
