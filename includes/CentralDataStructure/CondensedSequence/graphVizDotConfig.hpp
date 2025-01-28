#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP

#include "includes/CodeUtils/directories.hpp"

#include <filesystem>

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

} // namespace cdsCondensedSequence
#endif
