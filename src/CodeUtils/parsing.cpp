#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <string>
#include <algorithm>
#include <optional>
#include <stdexcept>

namespace codeUtils
{
    std::optional<int> parseInt(const std::string& str)
    {
        bool allDigits           = std::all_of(str.begin(), str.end(), ::isdigit);
        bool firstDashRestDigits = startsWith(str, "-") && std::all_of(str.begin() + 1, str.end(), ::isdigit);
        if (!(allDigits || firstDashRestDigits))
        {
            return std::nullopt;
        }
        try
        {
            return std::stoi(str);
        }
        catch (...)
        {
            return std::nullopt;
        }
        return int(0);
    }

    std::optional<ulong> parseUlong(const std::string& str)
    {
        if (!std::all_of(str.begin(), str.end(), ::isdigit))
        {
            return std::nullopt;
        }
        try
        {
            return std::stoul(str);
        }
        catch (...)
        {
            return std::nullopt;
        }
        return ulong(0);
    }

    std::optional<bool> parseBool(const std::string& str)
    {
        if (str == "true")
        {
            return true;
        }
        else if (str == "false")
        {
            return false;
        }
        else
        {
            return std::nullopt;
        }
        return false;
    };
} // namespace codeUtils
