#include "includes/CodeUtils/containers.hpp"

std::string codeUtils::FindStringInStringMap(const std::string& s,
                                             const std::unordered_map<std::string, std::string>& sMap)
{
    auto found = sMap.find(s);
    return (found != sMap.end()) ? found->second : "";
}

std::vector<bool> codeUtils::vectorAnd(const std::vector<bool>& vecA, const std::vector<bool>& vecB)
{
    std::vector<bool> result(vecA.size(), false);
    for (size_t n = 0; n < vecA.size(); n++)
    {
        result[n] = vecA[n] && vecB[n];
    }
    return result;
}

std::vector<size_t> codeUtils::offsetIndices(size_t offset, std::vector<size_t> indices)
{
    for (auto& index : indices)
    {
        index += offset;
    }
    return indices;
}

std::vector<bool> codeUtils::indicesToBools(size_t size, const std::vector<size_t>& indices)
{
    std::vector<bool> result(size, false);
    for (size_t n : indices)
    {
        result[n] = true;
    }
    return result;
}

std::vector<size_t> codeUtils::boolsToIndices(const std::vector<bool>& mask)
{
    size_t count = 0;
    for (bool b : mask)
    {
        count += b;
    }
    std::vector<size_t> result;
    result.reserve(count);
    for (size_t n = 0; n < mask.size(); n++)
    {
        if (mask[n])
        {
            result.push_back(n);
        }
    }
    return result;
}

std::vector<size_t> codeUtils::indexVectorWithOffset(size_t offset, size_t count)
{
    std::vector<size_t> result;
    result.reserve(count);
    for (size_t n = 0; n < count; n++)
    {
        result.push_back(n + offset);
    }
    return result;
}

std::vector<size_t> codeUtils::indexVector(size_t count)
{
    return indexVectorWithOffset(0, count);
}
