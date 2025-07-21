#include "include/util/parsing.hpp"

#include "include/util/strings.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>

namespace gmml
{
    namespace util
    {
        std::optional<int> parseInt(const std::string& str)
        {
            bool allDigits = std::all_of(str.begin(), str.end(), ::isdigit);
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

        std::optional<uint> parseUint(const std::string& str)
        {
            std::optional<int> parsed = parseInt(str);
            if (parsed.has_value())
            {
                int a = parsed.value();
                return (a < 0) ? std::nullopt : std::optional<uint> {uint(a)};
            }
            else
            {
                return std::nullopt;
            }
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

        std::optional<double> parseDouble(const std::string& str)
        {
            try
            {
                return std::stod(str);
            }
            catch (...)
            {
                return std::nullopt;
            }
            return double(0);
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
    } // namespace util
} // namespace gmml
