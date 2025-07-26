#ifndef INCLUDE_UTIL_STRINGS_HPP
#define INCLUDE_UTIL_STRINGS_HPP

#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace gmml
{
    namespace util
    {
        struct Brackets
        {
            std::string open;
            std::string close;
        };

        static const Brackets noBrackets = {"", ""};
        static const Brackets squareBrackets = {"[", "]"};
        static const Brackets roundBrackets = {"(", ")"};
        static const Brackets curlyBrackets = {"{", "}"};

        bool startsWith(const std::string& bigString, const std::string& smallString);
        std::string RemoveWhiteSpace(std::string str);
        void RemoveSpaces(std::string& str);
        void RemoveQuotes(std::string& str);
        std::string withoutSpaces(std::string str);
        std::string withoutQuotes(std::string str);
        int GetSizeOfIntInString(const std::string& str);
        std::string& Trim(std::string& str);
        void removeMultipleSpaces(std::string& str);
        void trimLeft(std::function<bool(const char&)> condition, std::string& str);
        void trimRight(std::function<bool(const char&)> condition, std::string& str);
        void trimWhitespace(std::string& str);
        std::string trimmedOfWhitespace(std::string str);
        std::vector<std::string> split(const std::string& s, char delim);
        std::string join(const std::string& delim, const std::vector<std::string>& strings);
        std::string upperCase(std::string str);
        std::string lowerCase(std::string str);
        std::string truncate(size_t n, const std::string& str);

        template<typename T> T from_string(const std::string& str)
        {
            std::istringstream iss(str);
            T value;
            iss >> value;
            return value;
        }

        template<class T> inline T extractFromStream(std::istream& ss, T temp)
        {
            ss >> temp;
            return temp;
        }
    } // namespace util
} // namespace gmml
#endif
