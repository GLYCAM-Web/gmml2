#include "include/util/filesystem.hpp"

#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace gmml
{
    namespace util
    {
        Path toPath(const std::string& str)
        {
            bool absolute = str[0] == '/';
            return Path {absolute, split(str, '/')};
        }

        std::string toString(const Path& path) { return (path.absolute ? "/" : "") + join("/", path.constituents); }

        bool pathExists(const std::string& path) { return std::filesystem::exists(path); }

        bool isDirectory(const std::string& path)
        {
            struct stat info;
            stat(path.c_str(), &info);
            return info.st_mode & S_IFDIR;
        }

        void ensureFileExists(const std::string& fileName)
        {
            if (!pathExists(fileName))
            {
                log(__LINE__, __FILE__, ERR, "File " + fileName + " does not exist");
                throw std::runtime_error("File " + fileName + " does not exist");
            }
        }

        bool directoryExists(const std::string& path) { return pathExists(path) && isDirectory(path); }

        bool directoryIsEmptyOrNonexistent(const std::string& pathName)
        {
            return !directoryExists(pathName) || std::filesystem::is_empty(pathName);
        }

        void createDirectories(const std::string& pathName)
        {
            std::filesystem::create_directories(pathName);
            struct stat info;
            if (stat(pathName.c_str(), &info) != 0)
            {
                throw std::runtime_error("failed to create directory: " + pathName);
            }
        }

        std::string readSymlink(const std::string& filename)
        {
            return std::filesystem::is_symlink(filename) ? std::filesystem::read_symlink(filename) : "";
        }

        Path pathToCurrentExecutable()
        {
            Path path = toPath(readSymlink("/proc/" + std::to_string(getpid()) + "/exe"));
            log(__LINE__, __FILE__, INF, "program running from: " + toString(path));
            return path;
        }

        Path pathAboveCurrentExecutableDir()
        {
            Path path = pathToCurrentExecutable();
            size_t pathSize = path.constituents.size();
            if (pathSize < 3)
            {
                throw std::runtime_error("Unexpected application path: " + join("/", path.constituents));
            }
            return {path.absolute, take(pathSize - 2, path.constituents)};
        }

        std::string SNFGSymbolsDir() { return "SNFG"; }
    } // namespace util
} // namespace gmml
