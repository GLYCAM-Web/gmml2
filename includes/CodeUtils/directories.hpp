#ifndef INCLUDES_CODEUTILS_DIRECTORIES_HPP
#define INCLUDES_CODEUTILS_DIRECTORIES_HPP
#include <sys/stat.h>  // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include <string>
#include <vector>
#include <filesystem>

// #include "includes/CodeUtils/logging.hpp"

namespace codeUtils
{

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

    struct Path
    {
        bool absolute;
        std::vector<std::string> constituents;
    };

    Path toPath(const std::string& str);
    bool doesDirectoryExist(const std::string& pathName);
    bool directoryIsEmptyOrNonexistent(const std::string& pathName);
    void createDirectories(const std::string& pathName);
    std::string getGmmlHomeDir();
    std::string getSNFGSymbolsDir();

    // TODO: Fix this(ese) abomination(s)
    static const std::filesystem::path gmmlHomeDirPath =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().string() + "/";
    // NOTE: Gotta double up the parent path to get rid of the
    // very last slash....
    static const std::filesystem::path gemsHomeDirPath =
        std::filesystem::path(gmmlHomeDirPath).parent_path().parent_path().string() + "/";

    static const std::filesystem::path relativeSnfgSymbolDirPath =
        "includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/";

    static const std::filesystem::path snfgSymbolsDirPath =
        std::filesystem::path(gmmlHomeDirPath).string() + relativeSnfgSymbolDirPath.string();

} // namespace codeUtils
#endif
