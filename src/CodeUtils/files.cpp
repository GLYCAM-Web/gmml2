#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <istream>
#include <fstream>
#include <functional>
#include <stdexcept>

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

std::istream& codeUtils::safeGetline(std::istream& in, std::string& out)
{
    out.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.

    std::istream::sentry se(in, true);
    std::streambuf* sb = in.rdbuf();

    while (true)
    {
        int c = sb->sbumpc();
        switch (c)
        {
            case '\n':
                return in;
            case '\r':
                if (sb->sgetc() == '\n')
                {
                    sb->sbumpc();
                }
                return in;
            case std::streambuf::traits_type::eof():
                // Also handle the case when the last line has no line ending
                if (out.empty())
                {
                    in.setstate(std::ios::eofbit);
                }
                return in;
            default:
                out += (char)c;
        }
    }
    return in;
}

void codeUtils::readFileLineByLine(const std::string& filename,
                                   std::function<void(const std::string&, size_t)> processLine)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::string message = filename + " could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    size_t lineNumber = 1;
    std::string out;
    while (!safeGetline(in, out).eof())
    {
        processLine(out, lineNumber);
        lineNumber++;
    }
}

std::vector<char> codeUtils::readEntireFile(const std::string& filename)
{
    ensureFileExists(filename);
    std::ifstream file(filename, std::ios::ate);
    // this works on text file with \n line endings on linux
    // might fail on other platforms, but the intent is to only use files we've generated ourselves
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size))
    {
        throw std::runtime_error("failed to read file " + filename);
    }
    return buffer;
}

void codeUtils::writeToFile(const std::string& filename, std::function<void(std::ostream&)> write)
{
    std::ofstream file(filename);
    write(file);
    file.close();
}
