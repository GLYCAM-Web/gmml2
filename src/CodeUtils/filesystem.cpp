#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/files.hpp"
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

std::string codeUtils::getSNFGSymbolsDir()
{
    return snfgSymbolsDirPath.string();
}

std::string codeUtils::getGmmlHomeDir()
{
    std::string gmmlHome = gmmlHomeDirPath.string();
    // TODO: Fix this gross ish
    if (gmmlHome.empty())
    {
        if (gemsHomeDirPath.string().empty())
        {
            std::string errorMessage =
                "$GMMLHOME and $GEMSHOME environmental variable not set (or std::getenv doesn't work on this system)";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
        gmmlHome = gemsHomeDirPath.string() + "/gmml/"; // guessing.
        if (!codeUtils::directoryExists(gmmlHome))
        {
            std::string errorMessage = "$GMMLHOME not set and directory $GEMSHOME/gmml/ doesn't exist.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
    }
    return gmmlHome;
}
