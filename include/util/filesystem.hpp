#ifndef INCLUDES_CODEUTILS_FILESYSTEM_HPP
#define INCLUDES_CODEUTILS_FILESYSTEM_HPP

#include <filesystem>
#include <string>
#include <vector>

namespace gmml
{
    namespace util
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
    } // namespace util
} // namespace gmml
#endif
