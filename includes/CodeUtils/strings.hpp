#ifndef INCLUDES_CODEUTILS_STRINGS_HPP
#define INCLUDES_CODEUTILS_STRINGS_HPP

#include <string>
#include <vector>
#include <sstream>

namespace codeUtils
{
    // Functions
    bool startsWith(std::string bigString, std::string smallString);
    std::string RemoveWhiteSpace(std::string str);
    void RemoveSpaces(std::string& str);
    void RemoveQuotes(std::string& str);
    int GetSizeOfIntInString(const std::string str);
    std::string& Trim(std::string& str);
    void removeMultipleSpaces(std::string& str);
    std::vector<std::string> split(const std::string& s, char delim);
    std::string join(const std::string& delim, const std::vector<std::string>& strings);

    template<typename T> T from_string(const std::string& str)
    {
        std::istringstream iss(str);
        T value;
        iss >> value;
        return value;
    }

    template<typename Out> inline void split(const std::string& s, char delim, Out result)
    {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim))
        {
            *(result++) = item;
        }
    }

    template<class T> inline T extractFromStream(std::istream& ss, T temp)
    {
        ss >> temp;
        return temp;
    }

} // namespace codeUtils
#endif
