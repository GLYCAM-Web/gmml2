#ifndef INCLUDES_CODEUTILS_FILES_HPP
#define INCLUDES_CODEUTILS_FILES_HPP

#include <string>

namespace codeUtils
{
    bool doesFileExist(const std::string& fileName);
    void ensureFileExists(const std::string& fileName);
    std::string SplitFilename(const std::string& str);
} // namespace codeUtils
#endif
