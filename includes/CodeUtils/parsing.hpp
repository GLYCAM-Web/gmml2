#ifndef INCLUDES_CODEUTILS_PARSING_HPP
#define INCLUDES_CODEUTILS_PARSING_HPP

#include <optional>
#include <string>

namespace codeUtils
{
    std::optional<int> parseInt(const std::string& str);
    std::optional<uint> parseUint(const std::string& str);
    std::optional<ulong> parseUlong(const std::string& str);
    std::optional<double> parseDouble(const std::string& str);
    std::optional<bool> parseBool(const std::string& str);
} // namespace codeUtils
#endif
