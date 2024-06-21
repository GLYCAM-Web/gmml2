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

    template<class T> void eraseNth(size_t n, std::vector<T>& vec)
    {
        vec.erase(vec.begin() + n);
    }

    template<class T> std::vector<T> withoutNth(size_t n, std::vector<T> vec)
    {
        eraseNth(n, vec);
        return vec;
    }

    template<class T> std::vector<T> vectorAppend(const std::vector<T> vecA, const std::vector<T> vecB)
    {
        std::vector<T> result;
        result.insert(result.end(), vecA.begin(), vecA.end());
        result.insert(result.end(), vecB.begin(), vecB.end());
        return result;
    }

    template<class T> std::vector<T*> pointerVector(std::vector<T>& vec)
    {
        std::vector<T*> result;
        result.reserve(vec.size());
        for (auto& a : vec)
        {
            result.push_back(&a);
        }
        return result;
    }

} // namespace codeUtils
#endif // GMML_INCLUDES_CODEUTILS_CONTAINERS_HPP
