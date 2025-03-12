#ifndef INCLUDES_CODEUTILS_FILES_HPP
#define INCLUDES_CODEUTILS_FILES_HPP

#include <cstddef>
#include <string>
#include <istream>
#include <functional>
#include <vector>

namespace codeUtils
{
    bool doesFileExist(const std::string& fileName);
    void ensureFileExists(const std::string& fileName);
    std::istream& safeGetline(std::istream& in, std::string& out);
    void readFileLineByLine(const std::string& filename, std::function<void(const std::string&, size_t)> processLine);
    std::vector<char> readEntireFile(const std::string& filename);
} // namespace codeUtils
#endif
