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

bool codeUtils::doesFileExist(const std::string& fileName)
{
    struct stat buffer;
    return (stat(fileName.c_str(), &buffer) == 0);
}

void codeUtils::ensureFileExists(const std::string& fileName)
{
    if (!codeUtils::doesFileExist(fileName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "File " + fileName + " does not exist");
        throw std::runtime_error("File " + fileName + " does not exist");
    }
}

bool codeUtils::doesDirectoryExist(const std::string& pathName)
{
    struct stat info;
    if (stat(pathName.c_str(), &info) != 0)
    {
        return false;
    }
    else if (info.st_mode & S_IFDIR)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool codeUtils::directoryIsEmptyOrNonexistent(const std::string& pathName)
{
    return !doesDirectoryExist(pathName) || std::filesystem::is_empty(pathName);
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
        if (!codeUtils::doesDirectoryExist(gmmlHome))
        {
            std::string errorMessage = "$GMMLHOME not set and directory $GEMSHOME/gmml/ doesn't exist.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
    }
    return gmmlHome;
}
