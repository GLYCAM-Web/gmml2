#include "include/util/strings.hpp"

#include <algorithm>
#include <ctype.h> // isdigit
#include <functional>
#include <sstream>

namespace gmml
{
    namespace util
    {
        bool startsWith(const std::string& bigString, const std::string& smallString)
        {
            return (bigString.compare(0, smallString.length(), smallString) == 0);
        }

        std::string RemoveWhiteSpace(std::string str)
        {
            RemoveSpaces(str);
            //    s.erase(s.find_last_not_of(" ") + 1);
            //    s.erase(0, s.find_first_not_of(" "));
            return str;
        }

        void RemoveSpaces(std::string& str) { str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end()); }

        void RemoveQuotes(std::string& str) { str.erase(std::remove(str.begin(), str.end(), '\"'), str.end()); }

        std::string withoutSpaces(std::string str)
        {
            RemoveSpaces(str);
            return str;
        }

        std::string withoutQuotes(std::string str)
        {
            RemoveQuotes(str);
            return str;
        }

        int GetSizeOfIntInString(const std::string& str)
        {
            int size = 0;
            for (const char& c : str)
            {
                if (isdigit(c))
                {
                    ++size;
                }
                else
                {
                    return size;
                }
            }
            return size;
        }

        std::string& Trim(std::string& str)
        {
            trimWhitespace(str);
            return str;
        }

        void removeMultipleSpaces(std::string& str)
        {
            size_t pos = str.find("  ");
            while (pos != std::string::npos)
            {
                str.erase(pos, 1);
                pos = str.find("  ");
            }
        }

        void trimLeft(std::function<bool(const char&)> condition, std::string& str)
        {
            str.erase(
                str.begin(),
                std::find_if(str.begin(), str.end(), [&condition](unsigned char ch) { return !condition(ch); }));
        }

        void trimRight(std::function<bool(const char&)> condition, std::string& str)
        {
            str.erase(
                std::find_if(str.rbegin(), str.rend(), [&condition](unsigned char ch) { return !condition(ch); })
                    .base(),
                str.end());
        }

        void trimWhitespace(std::string& str)
        {
            auto isSpace = [](const char c) { return std::isspace(c); };
            trimLeft(isSpace, str);
            trimRight(isSpace, str);
        }

        std::string trimmedOfWhitespace(std::string str)
        {
            trimWhitespace(str);
            return str;
        }

        std::vector<std::string> split(const std::string& s, char delim)
        {
            std::vector<std::string> result;
            std::stringstream ss(s);
            std::string item;
            while (std::getline(ss, item, delim))
            {
                if (item.size() >= 1 && !(item.size() == 1 && item[0] == delim))
                {
                    result.push_back(item);
                }
            }
            return result;
        }

        std::string join(const std::string& delim, const std::vector<std::string>& strings)
        {
            std::ostringstream result;
            for (size_t n = 0; n < strings.size(); n++)
            {
                if (n > 0)
                {
                    result << delim;
                }
                result << strings[n];
            }
            return result.str();
        }

        std::string upperCase(std::string str)
        {
            std::transform(str.begin(), str.end(), str.begin(), ::toupper);
            return str;
        }

        std::string lowerCase(std::string str)
        {
            std::transform(str.begin(), str.end(), str.begin(), ::tolower);
            return str;
        }

        std::string truncate(size_t n, const std::string& str) { return (str.size() <= n) ? str : str.substr(0, n); }
    } // namespace util
} // namespace gmml
