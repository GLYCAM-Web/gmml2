#include "includes/CodeUtils/find.hpp"

std::string codeUtils::FindStringInStringMap(const std::string& s,
                                             const std::unordered_map<std::string, std::string>& sMap)
{
    auto found = sMap.find(s);
    return (found != sMap.end()) ? found->second : "";
}
