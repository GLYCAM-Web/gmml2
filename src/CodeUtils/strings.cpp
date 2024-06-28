#include "includes/CodeUtils/strings.hpp"

#include <sstream>
#include <algorithm>
#include <iomanip> // std::setw
#include <ctype.h> // isdigit

bool codeUtils::startsWith(std::string bigString, std::string smallString)
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

void codeUtils::ExpandLine(std::string& line, int length)
{
    int l = line.length();
    if (l >= length)
    {
        return;
    }
    else
    {
        int space = length - l;
        std::stringstream ss;
        ss << line << std::setw(space) << " ";
        line = ss.str();
        return;
    }
    return;
}

int codeUtils::GetSizeOfIntInString(const std::string str)
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

#include <iostream>

std::string& codeUtils::Trim(std::string& str)
{
    if (str.size() > 2)
    {
        // Find last thing that isn't a space and delete from there to end.
        std::size_t lastNotSpacePosition = str.find_last_not_of(" ") + 1;
        if (lastNotSpacePosition < str.size())
        {
            str.erase(lastNotSpacePosition);
        }
        // Find first thing that isn't a space and delete from start to there.
        std::size_t firstNotSpacePosition = str.find_first_not_of(" ");
        if (firstNotSpacePosition > 0)
        {
            str.erase(0, str.find_first_not_of(" "));
        }
    }
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

std::vector<std::string> codeUtils::split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}
