#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <sys/param.h> // for MIN function
#include <sys/stat.h>  // stat
#include <cstring>     // strlen
#include <string>
#include <vector>

std::string codeUtils::Find_Program_Installation_Directory()
{ // A way to get the program name plus working directory
    char processID[32];
    char pBuffer[256];
    ssize_t len = sizeof(pBuffer);
    sprintf(processID, "/proc/%d/exe", getpid());
    int bytes = MIN(readlink(processID, pBuffer, len), len - 1);
    if (bytes >= 0)
    {
        pBuffer[bytes] = '\0';
        // std::cout << "processID:" << processID << " pBuffer:" << pBuffer << " bytes:" << bytes << std::endl;
        return codeUtils::SplitFilename(pBuffer);
    }
    return "Error";
}

std::string codeUtils::Find_Program_workingDirectory()
{
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
    {
        return "Error reading working directory";
    }
    cCurrentPath[strlen(cCurrentPath)] = '/';  // Add a / at the end.
    cCurrentPath[strlen(cCurrentPath)] = '\0'; // Above overwrites the null, the null is important. Respect the nu
    return cCurrentPath;
}

codeUtils::Path codeUtils::toPath(const std::string& str)
{
    bool absolute = str[0] == '/';
    return Path {absolute, split(str, '/')};
}

bool codeUtils::doesDirectoryExist(const std::string& pathName)
{
    struct stat info;
    if (stat(pathName.c_str(), &info) != 0)
    {
        return false; // printf( "cannot access %s\n", pathname );
    }
    else if (info.st_mode & S_IFDIR)
    {                // S_ISDIR() doesn't exist on my windows
        return true; // printf( "%s is a directory\n", pathname );
    }
    else
    {
        return false; // printf( "%s is no directory\n", pathname );
    }
}

void codeUtils::ensureDirectoryExists(const std::string& pathName)
{
    if (!codeUtils::doesDirectoryExist(pathName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Directory " + pathName + " does not exist");
        throw "Directory " + pathName + " does not exist";
    }
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

std::string codeUtils::getEnvVar(const std::string& key)
{
    char* val = std::getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
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
