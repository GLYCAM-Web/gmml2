#ifndef INCLUDES_CODEUTILS_FILES_HPP
#define INCLUDES_CODEUTILS_FILES_HPP

#include <string>
#include <istream>
#include <functional>

namespace codeUtils
{
    bool doesFileExist(const std::string& fileName);
    void ensureFileExists(const std::string& fileName);
    std::string SplitFilename(const std::string& str);
    std::istream& safeGetline(std::istream& in, std::string& out);
    void readFileLineByLine(std::string& filename, std::function<void(const std::string&)> processLine);
} // namespace codeUtils
#endif
