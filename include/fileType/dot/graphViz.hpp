#ifndef INCLUDE_FILETYPE_DOT_GRAPHVIZ_HPP
#define INCLUDE_FILETYPE_DOT_GRAPHVIZ_HPP

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace gmml
{
    namespace dot
    {
        struct Config
        {
            std::string basePath;
            std::string svgDirectory;
            std::string fileName;
            bool configLabels = true;
            bool edgeLabels = false;
            bool positionLabels = true;
            uint dpi = 72;
        };

        struct Image
        {
            bool found;
            std::string path;
            std::string label;
        };

        struct Edge
        {
            std::array<size_t, 2> nodeIndices;
            std::string label;
        };

        struct Node
        {
            size_t index;
            Image image;
            std::string label;
            std::string floatingLabel;
        };

        Image findImage(const Config& configs, const std::string& name, const std::string& label);
        std::string nodeBoxStyle(const Node& node);
        std::string nodeImageStyle(const Node& node);
        std::string edgeStyle(const std::vector<Node>& nodes, const Edge& linkage);
    } // namespace dot
} // namespace gmml

#endif
