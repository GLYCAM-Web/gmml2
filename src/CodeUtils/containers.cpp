#include "includes/CodeUtils/containers.hpp"

std::string codeUtils::FindStringInStringMap(const std::string& s,
                                             const std::unordered_map<std::string, std::string>& sMap)
{
    auto found = sMap.find(s);
    return (found != sMap.end()) ? found->second : "";
}

std::vector<bool> codeUtils::vectorAnd(const std::vector<bool>& vecA, const std::vector<bool>& vecB)
{
    std::vector<bool> result;
    result.reserve(vecA.size());
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
