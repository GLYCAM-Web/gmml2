#ifndef GMML_INCLUDES_CODEUTILS_CONTAINERS_HPP
#define GMML_INCLUDES_CODEUTILS_CONTAINERS_HPP

#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

namespace codeUtils
{
    std::string FindStringInStringMap(const std::string& s, const std::unordered_map<std::string, std::string>& sMap);

    template<class T> bool contains(const std::vector<T> vector, const T element)
    {
        return std::find(vector.begin(), vector.end(), element) != vector.end();
    }

} // namespace codeUtils
#endif // GMML_INCLUDES_CODEUTILS_CONTAINERS_HPP
