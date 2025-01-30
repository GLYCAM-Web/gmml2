#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZ_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZ_HPP

#include "includes/CodeUtils/directories.hpp"

#include <filesystem>
#include <string>
#include <array>
#include <vector>

namespace cdsCondensedSequence
{
    struct GraphVizDotConfig
    {
        GraphVizDotConfig(const std::string& basePath, const std::string& directoryPath, const std::string& filename)
            : file_name_(filename), base_path_(basePath), svg_directory_path_(directoryPath)
        {}

        bool show_config_labels_   = true;
        bool show_edge_labels_     = false;
        bool show_position_labels_ = true;
        int dpi_                   = 72;
        std::string file_name_;
        std::string base_path_;
        std::string svg_directory_path_;
    };

    struct GraphVizImage
    {
        bool found;
        std::string path;
        std::string label;
    };

    struct GraphVizLinkage
    {
        std::array<size_t, 2> nodeIndices;
        std::string label;
    };

    struct GraphVizResidueNode
    {
        size_t index;
        GraphVizImage image;
        std::string label;
        std::string floatingLabel;
    };

    GraphVizImage findImage(const cdsCondensedSequence::GraphVizDotConfig& configs, const std::string& name,
                            const std::string& label);
    std::string graphVizAglyconeNode(const GraphVizResidueNode& node);
    std::string graphVizSugarNode(const GraphVizResidueNode& node);
    std::string graphVizLinkageLine(const std::vector<GraphVizResidueNode>& nodes, const GraphVizLinkage& linkage);

} // namespace cdsCondensedSequence
#endif
