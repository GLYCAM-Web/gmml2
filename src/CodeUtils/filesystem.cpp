#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <sys/stat.h>
#include <string>
#include <vector>

codeUtils::Path codeUtils::toPath(const std::string& str)
{
    bool absolute = str[0] == '/';
    return Path {absolute, split(str, '/')};
}

std::string codeUtils::toString(const Path& path)
{
    return (path.absolute ? "/" : "") + codeUtils::join("/", path.constituents);
}

bool codeUtils::pathExists(const std::string& path)
{
    return std::filesystem::exists(path);
}

bool codeUtils::isDirectory(const std::string& path)
{
    struct stat info;
    stat(path.c_str(), &info);
    return info.st_mode & S_IFDIR;
}

void codeUtils::ensureFileExists(const std::string& fileName)
{
    if (!codeUtils::pathExists(fileName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "File " + fileName + " does not exist");
        throw std::runtime_error("File " + fileName + " does not exist");
    }
}

bool codeUtils::directoryExists(const std::string& path)
{
    return pathExists(path) && isDirectory(path);
}

bool codeUtils::directoryIsEmptyOrNonexistent(const std::string& pathName)
{
    return !directoryExists(pathName) || std::filesystem::is_empty(pathName);
}

void codeUtils::createDirectories(const std::string& pathName)
{
    std::filesystem::create_directories(pathName);
    struct stat info;
    if (stat(pathName.c_str(), &info) != 0)
    {
        throw std::runtime_error("failed to create directory: " + pathName);
    }
}

std::string codeUtils::readSymlink(const std::string& filename)
{
    return std::filesystem::is_symlink(filename) ? std::filesystem::read_symlink(filename) : "";
}

codeUtils::Path codeUtils::pathToCurrentExecutable()
{
    Path path = toPath(readSymlink("/proc/" + std::to_string(getpid()) + "/exe"));
    gmml::log(__LINE__, __FILE__, gmml::INF, "program running from: " + toString(path));
    return path;
}

codeUtils::Path codeUtils::pathAboveCurrentExecutableDir()
{
    Path path       = pathToCurrentExecutable();
    size_t pathSize = path.constituents.size();
    if (pathSize < 3)
    {
        throw std::runtime_error("Unexpected application path: " + join("/", path.constituents));
    }
    return {path.absolute, codeUtils::take(pathSize - 2, path.constituents)};
}

std::string codeUtils::SNFGSymbolsDir()
{
    return "includes/MolecularMetadata/Sugars/SNFG_Symbol_Images";
}
