#ifndef INCLUDES_CODEUTILS_CONTAINERS_HPP
#define INCLUDES_CODEUTILS_CONTAINERS_HPP

#include <algorithm>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

namespace codeUtils
{
    std::string FindStringInStringMap(const std::string& s, const std::unordered_map<std::string, std::string>& sMap);
    std::vector<bool> vectorAnd(const std::vector<bool>& vecA, const std::vector<bool>& vecB);
    std::vector<bool> vectorNot(const std::vector<bool>& vec);
    std::vector<size_t> offsetIndices(size_t offset, std::vector<size_t> indices);
    std::vector<bool> indicesToBools(size_t size, const std::vector<size_t>& indices);
    std::vector<size_t> boolsToIndices(const std::vector<bool>& mask);
    std::vector<bool> groupContains(const std::vector<size_t>& group, const std::vector<size_t>& indices);

    template<class T> size_t indexOf(const std::vector<T>& vector, const T element)
    {
        return std::find(vector.begin(), vector.end(), element) - vector.begin();
    }

    template<class T> std::vector<size_t> indicesOf(const std::vector<T>& vector, const std::vector<T>& elements)
    {
        std::vector<size_t> result;
        result.reserve(elements.size());
        for (auto& a : elements)
        {
            result.push_back(indexOf(vector, a));
        }
        return result;
    }

    template<class T> bool contains(const std::vector<T>& vector, const T element)
    {
        return indexOf(vector, element) < vector.size();
    }

    template<class T> std::vector<T> take(size_t count, const std::vector<T>& vec)
    {
        std::vector<T> result;
        count = std::min(count, vec.size());
        result.reserve(count);
        for (size_t n = 0; n < count; n++)
        {
            result.push_back(vec[n]);
        }
        return result;
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

    template<class T> std::vector<T> reverse(std::vector<T> vec)
    {
        std::reverse(vec.begin(), vec.end());
        return vec;
    }

    template<class A, class B> std::vector<A> vectorMap(std::function<A(const B&)>& func, const std::vector<B>& vec)
    {
        std::vector<A> result;
        result.reserve(vec.size());
        for (auto& a : vec)
        {
            result.push_back(func(a));
        }
        return result;
    }

    template<class T> std::vector<T> filter(std::function<bool(const T&)>& condition, const std::vector<T>& vec)
    {
        std::vector<T> result;
        result.reserve(vec.size());
        std::copy_if(vec.begin(), vec.end(), std::back_inserter(result), condition);
        return result;
    }

    template<class A, class T> A vectorSum(const A& initial, const std::vector<T>& vec)
    {
        A result = initial;
        for (size_t n = 0; n < vec.size(); n++)
        {
            result += vec[n];
        }
        return result;
    }

    template<class T> void insertInto(std::vector<T>& into, const std::vector<T>& other)
    {
        into.insert(into.end(), other.begin(), other.end());
    }

    template<class T> void fill(std::vector<T>& into, const T& value)
    {
        std::fill(into.begin(), into.end(), value);
    }

    template<class T> std::vector<T> vectorAppend(const std::vector<T>& vecA, const std::vector<T>& vecB)
    {
        std::vector<T> result;
        result.reserve(vecA.size() + vecB.size());
        insertInto(result, vecA);
        insertInto(result, vecB);
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

    template<class T> std::vector<T*> pointerToUniqueVector(std::vector<std::unique_ptr<T>>& vec)
    {
        std::vector<T*> result;
        result.reserve(vec.size());
        for (auto& a : vec)
        {
            result.push_back(a.get());
        }
        return result;
    }

    std::vector<size_t> indexVectorWithOffset(size_t offset, size_t count);
    std::vector<size_t> indexVector(size_t count);

    template<class T> std::vector<size_t> indexVectorWithOffset(size_t offset, const std::vector<T>& vec)
    {
        std::vector<size_t> result;
        result.reserve(vec.size());
        for (size_t n = 0; n < vec.size(); n++)
        {
            result.push_back(n + offset);
        }
        return result;
    }

    template<class T> std::vector<size_t> indexVector(const std::vector<T>& vec)
    {
        return indexVectorWithOffset(0, vec);
    }

    template<class T, class K>
    std::vector<T> indicesToValues(const std::vector<T>& values, const std::vector<K>& indices)
    {
        std::vector<T> result;
        result.reserve(indices.size());
        for (K n : indices)
        {
            result.push_back(values[n]);
        }
        return result;
    }

    template<class T> std::vector<T> boolsToValues(const std::vector<T>& values, const std::vector<bool>& mask)
    {
        std::vector<T> result;
        result.reserve(mask.size());
        for (size_t n = 0; n < mask.size(); n++)
        {
            if (mask[n])
            {
                result.push_back(values[n]);
            }
        }
        return result;
    }

    template<class T> std::vector<T> uniqueOnly(const std::vector<T>& values)
    {
        std::vector<T> result = values;
        std::sort(result.begin(), result.end());
        auto last = std::unique(result.begin(), result.end());
        result.erase(last, result.end());
        return result;
    }

    template<class T> std::vector<T> sorted(std::vector<T> values)
    {
        std::sort(values.begin(), values.end());
        return values;
    }

} // namespace codeUtils
#endif
