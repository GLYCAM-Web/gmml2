#ifndef INCLUDES_CODEUTILS_FILESYSTEM_HPP
#define INCLUDES_CODEUTILS_FILESYSTEM_HPP
#include <sys/stat.h>  // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include <unistd.h>
#include <string>
#include <vector>
#include <filesystem>

namespace codeUtils
{
    struct Path
    {
        bool absolute;
        std::vector<std::string> constituents;
    };

    Path toPath(const std::string& str);
    std::string toString(const Path& path);
    bool pathExists(const std::string& path);
    bool isDirectory(const std::string& path);
    void ensureFileExists(const std::string& fileName);
    bool directoryExists(const std::string& path);
    bool directoryIsEmptyOrNonexistent(const std::string& pathName);
    void createDirectories(const std::string& pathName);
    std::string readSymlink(const std::string& filename);
    Path pathToCurrentExecutable();
    Path pathAboveCurrentExecutableDir();
    std::string SNFGSymbolsDir();
    static const std::filesystem::path gmmlHomeDirPath =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().string();
} // namespace codeUtils
#endif
