#include "includes/CodeUtils/strings.hpp"

#include <sstream>
#include <algorithm>
#include <functional>
#include <ctype.h> // isdigit

bool codeUtils::startsWith(const std::string& bigString, const std::string& smallString)
{
    return (bigString.compare(0, smallString.length(), smallString) == 0);
}

std::string codeUtils::RemoveWhiteSpace(std::string str)
{
    codeUtils::RemoveSpaces(str);
    //    s.erase(s.find_last_not_of(" ") + 1);
    //    s.erase(0, s.find_first_not_of(" "));
    return str;
}

void codeUtils::RemoveSpaces(std::string& str)
{
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}

void codeUtils::RemoveQuotes(std::string& str)
{
    str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
}

std::string codeUtils::withoutSpaces(std::string str)
{
    RemoveSpaces(str);
    return str;
}

std::string codeUtils::withoutQuotes(std::string str)
{
    RemoveQuotes(str);
    return str;
}

int codeUtils::GetSizeOfIntInString(const std::string& str)
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

std::string& codeUtils::Trim(std::string& str)
{
    trimWhitespace(str);
    return str;
}

void codeUtils::removeMultipleSpaces(std::string& str)
{
    size_t pos = str.find("  ");
    while (pos != std::string::npos)
    {
        str.erase(pos, 1);
        pos = str.find("  ");
    }
}

void codeUtils::trimLeft(std::function<bool(const char&)> condition, std::string& str)
{
    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
                                        [&condition](unsigned char ch)
                                        {
                                            return !condition(ch);
                                        }));
}

void codeUtils::trimRight(std::function<bool(const char&)> condition, std::string& str)
{
    str.erase(std::find_if(str.rbegin(), str.rend(),
                           [&condition](unsigned char ch)
                           {
                               return !condition(ch);
                           })
                  .base(),
              str.end());
}

void codeUtils::trimWhitespace(std::string& str)
{
    auto isSpace = [](const char c)
    {
        return std::isspace(c);
    };
    trimLeft(isSpace, str);
    trimRight(isSpace, str);
}

std::string codeUtils::trimmedOfWhitespace(std::string str)
{
    trimWhitespace(str);
    return str;
}

std::vector<std::string> codeUtils::split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

std::string codeUtils::join(const std::string& delim, const std::vector<std::string>& strings)
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

std::string codeUtils::upperCase(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

std::string codeUtils::lowerCase(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

std::string codeUtils::truncate(size_t n, const std::string& str)
{
    n = std::min(n, str.size());
    return str.substr(0, n);
}
